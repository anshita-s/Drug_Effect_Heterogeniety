#!/usr/bin/env python3
"""
analyze_calibration.py
======================
Offline analysis of the files written by degradation_logger.py.  Runs OUTSIDE
CC3D with ordinary Python (needs numpy + matplotlib only).

Usage
-----
    python analyze_calibration.py RUNS_DIR [--ref-var 1.5] [--out analysis_out]
                                  [--checkpoints 100,200,300] [--bin 5]
                                  [--w-spatial 0.6 --w-temporal 0.2 --w-radial 0.1 --w-onset 0.1]
                                  [--w-coi 0.5]

RUNS_DIR is searched recursively for *_meta.json (one per run).

Comparison logic
----------------
* Reference runs = code_version "old" with var == --ref-var.
* Candidate runs = code_version "new".
* A candidate is only compared with references that have the SAME arrangement_id
  and the SAME matrix_signature (i.e. identical initial matrix cell ids/footprints).
  If a reference with the same seed exists it is used (paired); otherwise the
  candidate is compared with every reference of that arrangement and averaged.
* Spatial unit = ORIGINAL matrix cell ids.  "Degraded set at MCS t" = baseline
  cells that had a degradation event at or before t.  Because the layout code
  creates matrix cells in the same order in both versions, ids are comparable
  between old and new runs.  (Run id equality is checked through the baseline.)

Per component (COI, laminin) the score (lower = more similar) is
    w_spatial  * (1 - mean IoU over checkpoints)
  + w_temporal * min(1, RMSE of degraded-% curves / 100)
  + w_radial   * mean |difference of degraded fraction per radial bin| at the last MCS
  + w_onset    * min(1, |first-degradation MCS difference| / last MCS)
and the overall score is  w_coi*score_COI + (1-w_coi)*score_laminin.
Defaults put 60 % of the weight on spatial overlap.
"""
import os
import sys
import csv
import json
import glob
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMPS = ("COI", "LAMININ")
COLNAME = {"COI": "COI", "LAMININ": "laminin"}


# ----------------------------------------------------------------------------- loading
def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


class Run(object):
    def __init__(self, prefix):
        self.prefix = prefix
        with open(prefix + "_meta.json") as fh:
            self.meta = json.load(fh)
        self.code = self.meta["code_version"]
        self.var = float(self.meta["var"])
        self.arr = str(self.meta["arrangement_id"])
        self.seed = str(self.meta["seed"])
        self.sig = self.meta["matrix_signature"]
        self.final_mcs = self.meta.get("final_mcs")
        self.rows = read_csv(prefix + "_rows.csv")
        if self.final_mcs is None and self.rows:       # run did not close cleanly
            self.final_mcs = int(float(self.rows[-1]["MCS"]))
        events = read_csv(prefix + "_events.csv")
        base = read_csv(prefix + "_baseline_cells.csv")
        self.base = {}
        for c in COMPS:
            b = [r for r in base if r["component"] == c]
            self.base[c] = {
                "ids": np.array([int(r["cell_id"]) for r in b]),
                "xy": np.array([[fnum(r["x0"]), fnum(r["y0"])] for r in b]).reshape(-1, 2)}
        self.ev = {}
        for c in COMPS:
            e = [r for r in events if r["component"] == c and r["in_baseline"] == "1"]
            self.ev[c] = {"mcs": np.array([int(r["MCS"]) for r in e], dtype=int),
                          "id": np.array([int(r["cell_id"]) for r in e]),
                          "mmp": np.array([fnum(r["MMP_at_event"]) for r in e])}
        self.centroid = self.meta.get("cancer_centroid_baseline_xy") or [50.0, 50.0]

    def label(self):
        return "%s var=%g arr=%s seed=%s" % (self.code, self.var, self.arr, self.seed)

    def n0(self, c):
        return len(self.base[c]["ids"])

    def degraded_ids(self, c, t):
        e = self.ev[c]
        return set(e["id"][e["mcs"] <= t].tolist())

    def percent_curve(self, c, tgrid):
        m = np.sort(self.ev[c]["mcs"])
        n0 = self.n0(c)
        return 100.0 * np.searchsorted(m, tgrid, side="right") / n0 if n0 else np.zeros(len(tgrid))

    def first_mcs(self, c):
        m = self.ev[c]["mcs"]
        return int(m.min()) if len(m) else None


# ----------------------------------------------------------------------------- metrics
def radial_fraction(run, c, t, centroid, bin_w, rmax):
    xy = run.base[c]["xy"]
    r = np.hypot(xy[:, 0] - centroid[0], xy[:, 1] - centroid[1])
    deg = np.isin(run.base[c]["ids"], list(run.degraded_ids(c, t)))
    edges = np.arange(0, rmax + bin_w, bin_w)
    frac = np.full(len(edges) - 1, np.nan)
    for i in range(len(edges) - 1):
        sel = (r >= edges[i]) & (r < edges[i + 1])
        if sel.any():
            frac[i] = deg[sel].mean()
    return edges, frac


def compare(ref, cand, checkpoints, tend, args):
    out = {}
    centroid = ref.centroid
    for c in COMPS:
        if ref.n0(c) == 0:
            continue
        ious, precs, recs = [], [], []
        for t in checkpoints:
            a, b = ref.degraded_ids(c, t), cand.degraded_ids(c, t)
            u = len(a | b)
            if u == 0:
                continue                                   # nothing degraded in either run yet
            inter = len(a & b)
            ious.append(inter / u)
            recs.append(inter / len(a) if a else np.nan)   # share of reference-degraded cells also degraded in cand.
            precs.append(inter / len(b) if b else np.nan)  # share of cand.-degraded cells that are also in reference
        iou = float(np.mean(ious)) if ious else np.nan
        tg = np.arange(0, tend + 1)
        rmse = float(np.sqrt(np.mean((ref.percent_curve(c, tg) - cand.percent_curve(c, tg)) ** 2)))
        f_ref, f_c = ref.first_mcs(c), cand.first_mcs(c)
        if f_ref is None and f_c is None:
            onset = 0.0
        elif f_ref is None or f_c is None:
            onset = 1.0
        else:
            onset = min(1.0, abs(f_ref - f_c) / max(1, tend))
        _, fr = radial_fraction(ref, c, tend, centroid, args.bin, 80)
        _, fc = radial_fraction(cand, c, tend, centroid, args.bin, 80)
        ok = ~np.isnan(fr) & ~np.isnan(fc)
        radial = float(np.mean(np.abs(fr[ok] - fc[ok]))) if ok.any() else np.nan
        parts = {"spatial": 1.0 - iou if not np.isnan(iou) else 1.0,
                 "temporal": min(1.0, rmse / 100.0),
                 "radial": radial if not np.isnan(radial) else 1.0,
                 "onset": onset}
        score = (args.w_spatial * parts["spatial"] + args.w_temporal * parts["temporal"]
                 + args.w_radial * parts["radial"] + args.w_onset * parts["onset"])
        out[c] = {"iou": iou, "precision": float(np.nanmean(precs)) if len(precs) else np.nan,
                  "recall": float(np.nanmean(recs)) if len(recs) else np.nan,
                  "rmse_pct": rmse, "onset_diff": (None if f_ref is None or f_c is None else f_c - f_ref),
                  "final_pct_ref": float(ref.percent_curve(c, [tend])[0]),
                  "final_pct_cand": float(cand.percent_curve(c, [tend])[0]),
                  "radial_L1": radial, "score": score}
    sc = [out[c]["score"] for c in COMPS if c in out]
    if len(sc) == 2:
        out["overall"] = args.w_coi * out["COI"]["score"] + (1 - args.w_coi) * out["LAMININ"]["score"]
    elif sc:
        out["overall"] = sc[0]
    return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs_dir")
    ap.add_argument("--ref-var", type=float, default=1.5)
    ap.add_argument("--out", default="analysis_out")
    ap.add_argument("--checkpoints", default="")
    ap.add_argument("--bin", type=float, default=5.0, help="radial bin width (pixels)")
    ap.add_argument("--w-spatial", type=float, default=0.6)
    ap.add_argument("--w-temporal", type=float, default=0.2)
    ap.add_argument("--w-radial", type=float, default=0.1)
    ap.add_argument("--w-onset", type=float, default=0.1)
    ap.add_argument("--w-coi", type=float, default=0.5)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    metas = glob.glob(os.path.join(args.runs_dir, "**", "*_meta.json"), recursive=True)
    runs = []
    for m in sorted(metas):
        pre = m[:-len("_meta.json")]
        try:
            runs.append(Run(pre))
        except Exception as e:
            print("skipping %s (%s)" % (pre, e))
    if not runs:
        sys.exit("no runs found")
    refs = [r for r in runs if r.code == "old" and abs(r.var - args.ref_var) < 1e-12]
    cands = [r for r in runs if r.code == "new"]
    print("found %d runs: %d reference (old, var=%g), %d candidate (new)"
          % (len(runs), len(refs), args.ref_var, len(cands)))
    if not refs or not cands:
        sys.exit("need at least one reference and one candidate run")

    tend = int(min(r.final_mcs for r in runs))
    if args.checkpoints:
        checkpoints = [int(x) for x in args.checkpoints.split(",")]
    else:
        checkpoints = sorted(set(int(round(x)) for x in np.linspace(tend / 10.0, tend, 10)))
    print("common last MCS = %d; checkpoints = %s" % (tend, checkpoints))

    # suggest sweep range from the reference run(s)
    with open(os.path.join(args.out, "suggested_var_range.txt"), "w") as fh:
        for c in COMPS:
            mm = np.concatenate([r.ev[c]["mmp"] for r in refs]) if refs else np.array([])
            mm = mm[~np.isnan(mm)]
            if len(mm):
                q = np.percentile(mm, [5, 25, 50, 75, 95])
                line = "%s: MMP at old-code degradation events, percentiles 5/25/50/75/95 = %s" % (
                    c, ", ".join("%.4g" % v for v in q))
                print(line); fh.write(line + "\n")
        for r in sorted(runs, key=lambda r: (r.code, r.var, r.arr, r.seed)):
            mx = max([fnum(x["MMPc_max_COI"]) for x in r.rows if x["MMPc_max_COI"] != ""] or [np.nan])
            line = "max MMP at intact COI COM over run  [%s] = %.4g" % (r.label(), mx)
            fh.write(line + "\n")

    # compare every candidate with matching references
    per_run = []
    for cd in cands:
        pool = [r for r in refs if r.arr == cd.arr and r.sig == cd.sig]
        if not pool:
            print("WARNING: no reference with same arrangement/matrix for %s -- skipped" % cd.label())
            continue
        paired = [r for r in pool if r.seed == cd.seed]
        use = paired if paired else pool
        res = [compare(r, cd, checkpoints, tend, args) for r in use]
        row = {"arrangement": cd.arr, "var": cd.var, "seed": cd.seed,
               "paired_with_same_seed": int(bool(paired)), "n_refs": len(use)}
        for key in ("overall",):
            row[key] = float(np.mean([x[key] for x in res]))
        for c in COMPS:
            for k in ("iou", "precision", "recall", "rmse_pct", "final_pct_ref",
                      "final_pct_cand", "radial_L1", "score"):
                row["%s_%s" % (COLNAME[c], k)] = float(np.nanmean([x[c][k] for x in res]))
            fr, fc = use[0].first_mcs(c), cd.first_mcs(c)
            row["%s_first_MCS_ref" % COLNAME[c]] = fr
            row["%s_first_MCS_cand" % COLNAME[c]] = fc
        per_run.append(row)

    def write(name, rows):
        if not rows:
            return
        with open(os.path.join(args.out, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    write("summary_by_run.csv", per_run)

    # aggregate across seeds
    groups = {}
    for r in per_run:
        groups.setdefault((r["arrangement"], r["var"]), []).append(r)
    agg = []
    for (arr, v), rs in sorted(groups.items()):
        o = np.array([x["overall"] for x in rs])
        agg.append({"arrangement": arr, "var": v, "n_seeds": len(rs),
                    "overall_mean": o.mean(), "overall_sd": o.std(ddof=1) if len(o) > 1 else float("nan"),
                    "COI_iou_mean": np.nanmean([x["COI_iou"] for x in rs]),
                    "laminin_iou_mean": np.nanmean([x["laminin_iou"] for x in rs]),
                    "COI_rmse_pct": np.nanmean([x["COI_rmse_pct"] for x in rs]),
                    "laminin_rmse_pct": np.nanmean([x["laminin_rmse_pct"] for x in rs])})
    for arr in sorted(set(a["arrangement"] for a in agg)):
        sub = sorted([a for a in agg if a["arrangement"] == arr], key=lambda a: a["overall_mean"])
        rank = 1
        for a in sub:
            a["rank_in_arrangement"] = rank; rank += 1
    write("summary_by_var.csv", agg)

    # noise floor: reference seeds compared with each other
    floor = {}
    for arr in sorted(set(r.arr for r in refs)):
        rr = [r for r in refs if r.arr == arr]
        vals = []
        for i in range(len(rr)):
            for j in range(i + 1, len(rr)):
                if rr[i].sig == rr[j].sig:
                    vals.append(compare(rr[i], rr[j], checkpoints, tend, args)["overall"])
        if vals:
            floor[arr] = float(np.mean(vals))
    print("\nseed-to-seed noise floor of the score (old@%g vs old@%g, other seed): %s"
          % (args.ref_var, args.ref_var, floor if floor else "n/a (need >= 2 reference seeds)"))
    print("\n arrangement   var   seeds   overall(mean+-sd)   COI IoU  laminin IoU  rank")
    for a in sorted(agg, key=lambda a: (a["arrangement"], a["overall_mean"])):
        print(" %-11s %6g %6d   %8.3f +- %-6.3f   %7.3f   %9.3f   %3d" % (
            a["arrangement"], a["var"], a["n_seeds"], a["overall_mean"],
            0 if np.isnan(a["overall_sd"]) else a["overall_sd"],
            a["COI_iou_mean"], a["laminin_iou_mean"], a["rank_in_arrangement"]))

    # ------------------------------------------------------------------- plots
    for arr in sorted(set(a["arrangement"] for a in agg)):
        tg = np.arange(0, tend + 1)
        rr = [r for r in refs if r.arr == arr]
        for c in COMPS:
            fig, ax = plt.subplots(figsize=(7, 4.5))
            for v in sorted(set(r.var for r in cands if r.arr == arr)):
                cs = [r for r in cands if r.arr == arr and r.var == v]
                ax.plot(tg, np.mean([r.percent_curve(c, tg) for r in cs], axis=0), lw=1.4,
                        label="new var=%g (n=%d)" % (v, len(cs)))
            ax.plot(tg, np.mean([r.percent_curve(c, tg) for r in rr], axis=0), "k", lw=3,
                    label="OLD var=%g (n=%d)" % (args.ref_var, len(rr)))
            ax.set_xlabel("MCS"); ax.set_ylabel("%s degraded (%% of initial, cumulative)" % COLNAME[c])
            ax.set_title("Arrangement %s - %s" % (arr, COLNAME[c])); ax.grid(alpha=.3)
            ax.legend(fontsize=7)
            fig.tight_layout(); fig.savefig(os.path.join(args.out, "curve_%s_arr%s.png" % (COLNAME[c], arr)), dpi=150)
            plt.close(fig)

        sub = sorted([a for a in agg if a["arrangement"] == arr], key=lambda a: a["var"])
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.errorbar([a["var"] for a in sub], [a["overall_mean"] for a in sub],
                    yerr=[0 if np.isnan(a["overall_sd"]) else a["overall_sd"] for a in sub], marker="o")
        if arr in floor:
            ax.axhline(floor[arr], ls="--", c="gray", label="seed-to-seed floor (old vs old)")
            ax.legend(fontsize=8)
        ax.set_xscale("log"); ax.set_xlabel("new-code var"); ax.set_ylabel("overall score (lower = closer to old)")
        ax.set_title("Arrangement %s" % arr); ax.grid(alpha=.3)
        fig.tight_layout(); fig.savefig(os.path.join(args.out, "score_vs_var_arr%s.png" % arr), dpi=150)
        plt.close(fig)

        best = min(sub, key=lambda a: a["overall_mean"])
        cd = [r for r in cands if r.arr == arr and r.var == best["var"]][0]
        ref = ([r for r in rr if r.seed == cd.seed] or rr)[0]
        fig, axs = plt.subplots(1, 2, figsize=(11, 5))
        for ax, c in zip(axs, COMPS):
            xy = ref.base[c]["xy"]
            a = np.isin(ref.base[c]["ids"], list(ref.degraded_ids(c, tend)))
            b = np.isin(ref.base[c]["ids"], list(cd.degraded_ids(c, tend)))
            for mask, col, lab in ((~a & ~b, "#cccccc", "neither"), (a & b, "#2ca02c", "both"),
                                   (a & ~b, "#d62728", "old only"), (~a & b, "#1f77b4", "new only")):
                ax.scatter(xy[mask, 0], xy[mask, 1], s=14 if c == "COI" else 40, c=col, label=lab)
            ax.plot(*ref.centroid, "k+", ms=14)
            ax.set_aspect("equal"); ax.set_title("%s at MCS %d (arr %s, new var=%g)" % (COLNAME[c], tend, arr, best["var"]))
            ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(args.out, "overlay_best_arr%s.png" % arr), dpi=150)
        plt.close(fig)
    print("\nwrote results to %s/" % args.out)


if __name__ == "__main__":
    main()