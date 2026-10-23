"""
degradation_logger.py
=====================
Shared, model-neutral logger for matrix-degradation calibration in CompuCell3D.

Put this file in the SAME folder as your steppable file (the project's
Simulation/ folder).  Both the OLD (MMP/I) and NEW (MMP only) steppable files
import it, so the two versions use identical definitions, units and sampling.

The logger is READ-ONLY with respect to the model: it never changes cell types,
volumes, fields, cell dictionaries (CompuCell.getPyAttrib) or the random number
stream.  It only reads cell properties and writes files.

Definitions
-----------
* A "matrix cell" is a COI or LAMININ cell that exists (volume > 0) when the
  baseline is taken (right after CellLayoutSteppable.start()).
* A degradation EVENT is the moment MatrixDegradation switches a COI cell to
  C_LYSED, or a LAMININ cell to L_LYSED.  Each original matrix cell can have at
  most one event (a lysed cell never goes back to COI/LAMININ), and the logger
  also de-duplicates by cell id.
* Cumulative degraded count  = number of baseline cells that have had an event.
  It never decreases, and is not affected by the C_LYSED -> NCI conversion after
  40 MCS, nor by the empty "ghost" cells that MatrixDegradation.step() creates
  every MCS with potts.createCell().
* degraded_percent = 100 * cumulative degraded count / initial count.

Output files (all share one prefix, unique per run):
  <prefix>_rows.csv            one row per sampled MCS (+ final MCS)
  <prefix>_events.csv          one row per degradation event (spatial/temporal)
  <prefix>_baseline_cells.csv  every matrix cell at baseline: id, component, x0, y0, area0
  <prefix>_baseline_maps.npz   label_map (cell id per pixel), type_map at baseline
  <prefix>_types_mcsNNNNNN.npz optional snapshots of the type map (SNAPSHOT_EVERY)
  <prefix>_meta.json           run settings, signatures, initial counts, summary
"""
import os
import csv
import json
import time
import hashlib
import numpy as np


def _f(v):
    """Format a number for CSV; None -> empty field."""
    if v is None:
        return ""
    try:
        return "%.6g" % float(v)
    except (TypeError, ValueError):
        return str(v)


class DegradationLogger(object):
    COMPONENTS = ("COI", "LAMININ")
    _COL = {"COI": "COI", "LAMININ": "laminin"}   # column-name stem

    ROW_FIELDS = (
        ["MCS", "var", "code_version", "arrangement_id", "seed"]
        + sum([[
            "initial_%s_count" % s,
            "intact_%s_count" % s,
            "%s_degraded_count" % s,
            "%s_degraded_percent" % s,
            "%s_currently_lysed" % s,
            "%s_lost_undegraded" % s,
        ] for s in ("COI", "laminin")], [])
        + ["first_COI_degradation_MCS", "first_laminin_degradation_MCS",
           "MMPc_max_COI", "MMPc_max_laminin"]
    )

    EVENT_FIELDS = ["MCS", "component", "cell_id", "in_baseline",
                    "x_event", "y_event", "volume_at_event",
                    "MMP_at_event", "I_at_event", "T1_at_event",
                    "nearest_cancer_dist"]

    def __init__(self, code_version, var, arrangement_id, seed,
                 log_every=10, snapshot_every=100, out_dir=None):
        self.code_version = str(code_version)
        self.var = float(var)
        self.arrangement_id = str(arrangement_id)
        self.seed = seed
        self.log_every = max(1, int(log_every))
        self.snapshot_every = int(snapshot_every)
        self.out_dir_cfg = out_dir

        self.baseline_done = False
        self.baseline_mcs = None
        self.closed = False

        self.initial = dict((c, 0) for c in self.COMPONENTS)
        self.baseline_ids = dict((c, set()) for c in self.COMPONENTS)
        self.degraded = dict((c, {}) for c in self.COMPONENTS)   # id -> MCS
        self._seen = dict((c, set()) for c in self.COMPONENTS)
        self.extra_events = dict((c, 0) for c in self.COMPONENTS)
        self.first_mcs = dict((c, None) for c in self.COMPONENTS)
        self._mmp_max = dict((c, None) for c in self.COMPONENTS)

        self._cancer_xy = np.zeros((0, 2))
        self._last_mcs = None
        self._last_row_mcs = None
        self._last_snap_mcs = None
        self._warned_step = False
        self._dims = None

        self._row_fh = self._row_w = self._evt_fh = self._evt_w = None
        self.prefix = None
        self.meta = {}

    # ------------------------------------------------------------------ helpers
    def percent(self, comp):
        n0 = self.initial.get(comp, 0)
        return 100.0 * len(self.degraded[comp]) / n0 if n0 else 0.0

    @staticmethod
    def _get_dims(steppable):
        try:
            d = steppable.dim
            return int(d.x), int(d.y)
        except Exception:
            d = steppable.cell_field.getDim()
            return int(d.x), int(d.y)

    def _resolve_out_dir(self, steppable):
        env = os.environ.get("CALIB_OUTDIR")
        if env:
            d = env
        elif self.out_dir_cfg:
            d = self.out_dir_cfg
        else:
            base = getattr(steppable, "output_dir", None)
            if not base:
                base = os.path.join(os.getcwd(), "degradation_output")
            d = os.path.join(str(base), "degradation_logs")
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

    def _lattice_maps(self, steppable):
        nx, ny = self._dims
        label = np.zeros((nx, ny), dtype=np.int32)
        ctype = np.zeros((nx, ny), dtype=np.int16)
        cf = steppable.cell_field
        for x in range(nx):
            for y in range(ny):
                c = cf[x, y, 0]
                if c is not None:
                    label[x, y] = c.id
                    ctype[x, y] = c.type
        return label, ctype

    # ----------------------------------------------------------------- baseline
    def ensure_baseline(self, steppable, mcs):
        """Take the baseline once the matrix exists.  Safe to call repeatedly."""
        if self.baseline_done:
            return True
        coi_all = list(steppable.cell_list_by_type(steppable.COI))
        lam_all = list(steppable.cell_list_by_type(steppable.LAMININ))
        coi = [c for c in coi_all if c.volume > 0]
        lam = [c for c in lam_all if c.volume > 0]
        if not coi and not lam:
            return False                    # layout not created yet

        self._dims = self._get_dims(steppable)
        self.baseline_mcs = int(mcs)
        out_dir = self._resolve_out_dir(steppable)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        self.prefix = os.path.join(
            out_dir, "degr_%s_arr%s_var%g_seed%s_%s" % (
                self.code_version, self.arrangement_id, self.var, self.seed, stamp))

        base_rows = []
        for comp, cells in (("COI", coi), ("LAMININ", lam)):
            self.baseline_ids[comp] = set(c.id for c in cells)
            self.initial[comp] = len(cells)
            for c in cells:
                base_rows.append([c.id, comp, _f(c.xCOM), _f(c.yCOM), _f(c.volume)])
        with open(self.prefix + "_baseline_cells.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cell_id", "component", "x0", "y0", "area0"])
            w.writerows(base_rows)

        label, ctype = self._lattice_maps(steppable)
        all_ids = np.array(sorted(self.baseline_ids["COI"] | self.baseline_ids["LAMININ"]),
                           dtype=np.int32)
        matrix_label = np.where(np.isin(label, all_ids), label, 0).astype(np.int32)
        np.savez_compressed(self.prefix + "_baseline_maps.npz",
                            label_map=label, type_map=ctype)

        cancer = list(steppable.cell_list_by_type(steppable.CANCER_S, steppable.CANCER_R))
        if cancer:
            cx = float(np.mean([c.xCOM for c in cancer]))
            cy = float(np.mean([c.yCOM for c in cancer]))
        else:
            cx = cy = None

        self.meta = {
            "code_version": self.code_version, "var": self.var,
            "arrangement_id": self.arrangement_id, "seed": self.seed,
            "log_every": self.log_every, "snapshot_every": self.snapshot_every,
            "lattice_xy": list(self._dims), "baseline_mcs": self.baseline_mcs,
            "initial_COI_count": self.initial["COI"],
            "initial_laminin_count": self.initial["LAMININ"],
            "zero_volume_matrix_cells_excluded": (len(coi_all) - len(coi)) + (len(lam_all) - len(lam)),
            "type_ids": dict((n, int(getattr(steppable, n))) for n in
                             ("COI", "LAMININ", "C_LYSED", "L_LYSED", "NCI",
                              "CANCER_S", "CANCER_R")),
            "matrix_signature": hashlib.md5(matrix_label.tobytes()).hexdigest()[:12],
            "layout_signature": hashlib.md5(ctype.tobytes()).hexdigest()[:12],
            "cancer_centroid_baseline_xy": [cx, cy],
            "start_time": stamp,
            "note": "Potts RNG seed is set in the XML (<RandomSeed>), not by this logger.",
        }
        self._write_meta()

        self._row_fh = open(self.prefix + "_rows.csv", "w", newline="")
        self._row_w = csv.DictWriter(self._row_fh, fieldnames=self.ROW_FIELDS)
        self._row_w.writeheader()
        self._evt_fh = open(self.prefix + "_events.csv", "w", newline="")
        self._evt_w = csv.DictWriter(self._evt_fh, fieldnames=self.EVENT_FIELDS)
        self._evt_w.writeheader()
        self._row_fh.flush()
        self._evt_fh.flush()
        self.baseline_done = True
        print("[DegradationLogger] baseline at MCS %s: COI=%d LAMININ=%d -> %s_*"
              % (mcs, self.initial["COI"], self.initial["LAMININ"], self.prefix))
        return True

    def _write_meta(self):
        with open(self.prefix + "_meta.json", "w") as fh:
            json.dump(self.meta, fh, indent=2, sort_keys=True)

    # --------------------------------------------------------- per-step hooks
    def begin_step(self, steppable, mcs):
        """Call at the START of MatrixDegradation.step()."""
        self.ensure_baseline(steppable, mcs)
        if not self._warned_step:
            sim_step = steppable.simulator.getStep()
            if sim_step != mcs:
                print("[DegradationLogger] WARNING: step arg mcs=%s but simulator.getStep()=%s"
                      % (mcs, sim_step))
            self._warned_step = True
        cancer = list(steppable.cell_list_by_type(steppable.CANCER_S, steppable.CANCER_R))
        if cancer:
            self._cancer_xy = np.array([[c.xCOM, c.yCOM] for c in cancer], dtype=float)
        else:
            self._cancer_xy = np.zeros((0, 2))

    def observe_mmp(self, comp, value):
        """Track the largest MMP value seen at intact matrix cell COMs since the last row."""
        v = float(value)
        cur = self._mmp_max[comp]
        if cur is None or v > cur:
            self._mmp_max[comp] = v

    def record_event(self, mcs, comp, cell, mmp, i_val, t1):
        """Call right after a matrix cell has been switched to its lysed type."""
        if not self.baseline_done:
            return
        cid = cell.id
        if cid in self._seen[comp]:
            return                                   # never double count
        self._seen[comp].add(cid)
        in_base = cid in self.baseline_ids[comp]
        if in_base:
            self.degraded[comp][cid] = int(mcs)
        else:
            self.extra_events[comp] += 1
        if self.first_mcs[comp] is None:
            self.first_mcs[comp] = int(mcs)
        x, y = float(cell.xCOM), float(cell.yCOM)
        if len(self._cancer_xy):
            d = np.sqrt(((self._cancer_xy - np.array([x, y])) ** 2).sum(axis=1)).min()
        else:
            d = None
        self._evt_w.writerow({
            "MCS": int(mcs), "component": comp, "cell_id": cid,
            "in_baseline": int(in_base), "x_event": _f(x), "y_event": _f(y),
            "volume_at_event": _f(cell.volume), "MMP_at_event": _f(mmp),
            "I_at_event": _f(i_val), "T1_at_event": _f(t1),
            "nearest_cancer_dist": _f(d)})

    def end_step(self, steppable, mcs):
        """Call at the END of MatrixDegradation.step()."""
        if not self.baseline_done:
            return
        self._last_mcs = int(mcs)
        if mcs % self.log_every == 0:
            self.log_row(steppable, mcs)
        if self.snapshot_every > 0 and mcs % self.snapshot_every == 0:
            self.snapshot(steppable, mcs)

    # ------------------------------------------------------------- output
    def log_row(self, steppable, mcs):
        row = {"MCS": int(mcs), "var": self.var, "code_version": self.code_version,
               "arrangement_id": self.arrangement_id, "seed": self.seed}
        for comp, ctype, ltype in (("COI", steppable.COI, steppable.C_LYSED),
                                   ("LAMININ", steppable.LAMININ, steppable.L_LYSED)):
            s = self._COL[comp]
            ids = self.baseline_ids[comp]
            intact = sum(1 for c in steppable.cell_list_by_type(ctype) if c.id in ids)
            lysed_now = sum(1 for c in steppable.cell_list_by_type(ltype)
                            if c.id in self.degraded[comp])
            n0 = self.initial[comp]
            deg = len(self.degraded[comp])
            row["initial_%s_count" % s] = n0
            row["intact_%s_count" % s] = intact
            row["%s_degraded_count" % s] = deg
            row["%s_degraded_percent" % s] = _f(100.0 * deg / n0) if n0 else ""
            row["%s_currently_lysed" % s] = lysed_now
            row["%s_lost_undegraded" % s] = n0 - intact - deg
        row["first_COI_degradation_MCS"] = "" if self.first_mcs["COI"] is None else self.first_mcs["COI"]
        row["first_laminin_degradation_MCS"] = "" if self.first_mcs["LAMININ"] is None else self.first_mcs["LAMININ"]
        row["MMPc_max_COI"] = _f(self._mmp_max["COI"])
        row["MMPc_max_laminin"] = _f(self._mmp_max["LAMININ"])
        self._mmp_max = dict((c, None) for c in self.COMPONENTS)
        self._row_w.writerow(row)
        self._row_fh.flush()
        self._evt_fh.flush()
        self._last_row_mcs = int(mcs)

    def snapshot(self, steppable, mcs):
        if self._last_snap_mcs == mcs:
            return
        _, ctype = self._lattice_maps(steppable)
        np.savez_compressed("%s_types_mcs%06d.npz" % (self.prefix, mcs),
                            type_map=ctype, mcs=int(mcs))
        self._last_snap_mcs = int(mcs)

    def close(self, steppable):
        """Idempotent.  Called from finish() / on_stop()."""
        if self.closed or not self.baseline_done:
            return
        self.closed = True
        if self._last_mcs is not None:
            if self._last_row_mcs != self._last_mcs:
                self.log_row(steppable, self._last_mcs)          # guarantee a final row
            if self.snapshot_every > 0:
                self.snapshot(steppable, self._last_mcs)
        for fh in (self._row_fh, self._evt_fh):
            if fh is not None:
                fh.flush()
                fh.close()
        self.meta.update({
            "final_mcs": self._last_mcs,
            "first_COI_degradation_MCS": self.first_mcs["COI"],
            "first_laminin_degradation_MCS": self.first_mcs["LAMININ"],
            "final_COI_degraded_percent": self.percent("COI"),
            "final_laminin_degraded_percent": self.percent("LAMININ"),
            "events_not_in_baseline": dict(self.extra_events),
            "finished_cleanly": True})
        self._write_meta()
        print("[DegradationLogger] closed: COI %.1f%%, laminin %.1f%% degraded at MCS %s"
              % (self.percent("COI"), self.percent("LAMININ"), self._last_mcs))