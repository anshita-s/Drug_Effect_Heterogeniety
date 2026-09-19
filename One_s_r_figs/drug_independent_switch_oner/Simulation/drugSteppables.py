# Importing python libraries
from cc3d.cpp.PlayerPython import * 
from cc3d.core.PySteppables import *
from cc3d.cpp import CompuCell
import sys
from cc3d import CompuCellSetup
import numpy as np
import random
from random import uniform
import math
from math import *
import string
from os import getcwd
import os.path
from xml.dom import minidom
import itertools as it


#location = os.path.dirname(os.path.realpath(__file__))

# Few global variables
var = 1.5 #RD_coop value
varii_1 = 2.25 #growth rate coefficient 
varii_2 = 1

c = 0
mcsOut = 0

# --- Drug-response parameters ---
DRUG_DEATH_RATE = 1.0    # targetVolume lost per MCS once a CANCER_S cell is 'hit'
DRUG_MIN_VOLUME = 2.0    # once volume drops below this, the cell is removed

# --- Drug uptake parameters --
DRUG_UPTAKE_RATE_S  = 0.05   # fraction of local [DRUG] absorbed per MCS by CANCER_S
DRUG_UPTAKE_RATE_R  = 0.01   # fraction of local [DRUG] absorbed per MCS by CANCER_R
DRUG_UPTAKE_THRESHOLD = 1.5  # cumulative drugUptake that triggers the fate decision in CANCER_S

# --- Random S<->R type-switch parameters ---
# Independent of drug threshold and of mitosis. Every so often, a random
# (x, y) lattice site is picked and whatever CANCER_S/CANCER_R cell sits
# there has its type flipped. The wait between events is drawn from
# Normal(MU_SWITCH_WINDOW, SIGMA_SWITCH_WINDOW), redrawn after every event.
MU_SWITCH_WINDOW    = 100.0  # mean wait between switch events (MCS)
SIGMA_SWITCH_WINDOW = 15.0   # std-dev of wait between switch events (MCS)
MIN_SWITCH_WINDOW   = 1.0    # floor, so a sampled wait can never be <= 0

# --- Forced S+S division after an S->R switch ---
# Whenever RandomTypeSwitchSteppable converts a cell S -> R, this flag is
# raised. Whichever cell divides *next* anywhere in the lattice -- whether
# its own parent is CANCER_S or CANCER_R -- must produce two CANCER_S
# daughters. Both MitosisSteppable_S and MitosisSteppable_R consume this
# flag (whichever fires first), so it applies to exactly one division.
FORCE_NEXT_DIVISION_SS = False


def _consume_force_next_division_ss():
    """If a forced-S+S division is pending, clear it and return True."""
    global FORCE_NEXT_DIVISION_SS
    if FORCE_NEXT_DIVISION_SS:
        FORCE_NEXT_DIVISION_SS = False
        return True
    return False



############################# INITIAL CELLS AND ECM ARRANGEMENT #########################
class CellLayoutSteppable(SteppableBasePy):
    
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
    def start(self):  
        
        #COLLAGEN-I
        for x in range(5,95,4):
            for y in range(5,95,4):
                self.cell_field[x:x+2,y:y+2,0] = self.new_cell(self.COI)
                
        #LAMININ      
        for x in range(40,60,3):
            for y in range(40,60,3):
                self.cell_field[x:x+3,y:y+3,0] = self.new_cell(self.LAMININ)        
                
                
    #    # Drug-SENSITIVE subset (formerly CANCER_A, "Low HI")
    #    for x in range(42,58,4):
    #        for y in range(50,58,4):
    #            self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_S)
    #           
    #    # Drug-RESISTANT subset (formerly CANCER_B)
    #    for x in range(42,58,4):
    #        for y in range(42,50,4):
    #            self.cell_field[x:x+4, y:y+4,0] = self.new_cell(self.CANCER_R)
                
        
        # # # High HI
     #   for x in it.chain(range(42,46,4), range(50,54,4)):
     #       for y in it.chain(range(42,46,4), range(50,54,4)):
     #           self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_S)
     #   for x in it.chain(range(46,50,4), range(54,58,4)):
     #       for y in it.chain(range(46,50,4), range(54,58,4)):
     #           self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_S)
     # 
     #   for x in it.chain(range(46,50,4), range(54,58,4)):
     #       for y in it.chain(range(42,46,4), range(50,54,4)):
     #           self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_R)
     #   for x in it.chain(range(42,46,4), range(50,54,4)):
     #       for y in it.chain(range(46,50,4), range(54,58,4)):
     #          self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_R)
      
     #   # # # Single-sensitive setup: 1 CANCER_S (bottom-left of the 4x4 grid), 15 CANCER_R
     #   x_positions = [42, 46, 50, 54]
     #   y_positions = [42, 46, 50, 54]
     #   sensitive_pos = (min(x_positions), min(y_positions))  # bottom-left corner cell

     #   for x in x_positions:
     #       for y in y_positions:
     #           if (x, y) == sensitive_pos:
     #              self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_S)
     #           else:
     #              self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_R)
                    
        # # # Single-sensitive setup: 1 CANCER_R (bottom-left of the 4x4 grid), 15 CANCER_S
        x_positions = [42, 46, 50, 54]
        y_positions = [42, 46, 50, 54]
        sensitive_pos = (min(x_positions), min(y_positions))  # bottom-left corner cell

        for x in x_positions:
            for y in y_positions:
                if (x, y) == sensitive_pos:
                    self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_R)
                else:
                    self.cell_field[x:x+4, y:y+4, 0] = self.new_cell(self.CANCER_S)                
                

class VolumeParamSteppable(SteppableBasePy):
    
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)


    def start(self):
        global c
        global ck
        global c_value
        global mcsOut
        global varii
        c_value =[]
        
        
        for cell in self.cell_list_by_type(self.COI):
            cell.targetVolume = 4.0
            cell.lambdaVolume = 20.0           
        
        for cell in self.cell_list_by_type(self.LAMININ):
            cell.targetVolume = 9.0
            cell.lambdaVolume = 20.0
        
        # for cell in self.cell_list_by_type(self.CANCER_S, self.CANCER_R):
            # cell.targetVolume = 24.0
            # cell.lambdaVolume = 20.0

        for cell in self.cell_list_by_type(self.CANCER_S):
            cell.targetVolume = 16.0 + 8.0 
            cell.lambdaVolume = 20.0
            cell.targetSurface = 20.0
            cell.lambdaSurface = 2.0  

        for cell in self.cell_list_by_type(self.CANCER_R):
            cell.targetVolume = 16.0 + 8.0 
            cell.lambdaVolume = 20.0
            cell.targetSurface = 20.0
            cell.lambdaSurface = 2.0 

    
    def step(self,mcs):
        
        # field = self.field.FIELD_NAME
        # value = field[10, 10, 0]        
        global c
        global ck
        global mcsOut
        global c_value
        global varii
        

        GF_field = self.field.GF
        
        for cell in self.cellList:      
            if cell.type == self.COI and mcs<10:
                cell.targetVolume+=0.8
                cell.lambdaVolume=20.0      
        # for cell in self.cell_list_by_type(self.COI):
            # cell.targetVolume = 6.0
            # cell.lambdaVolume = 20.0      
     
            if cell.type == self.CANCER_S:
                cell.targetSurface = 20.0
                cell.lambdaSurface = 2.0
                
            if cell.type == self.CANCER_R:
                cell.targetSurface = 20.0
                cell.lambdaSurface = 2.0
        
        for cell in self.cell_list_by_type(self.CANCER_S):
            neighbor_list_1 = self.get_cell_neighbor_data_list(cell)
            k1 = neighbor_list_1.common_surface_area_with_cell_types(cell_type_list=[1,2])
            s1 = cell.surface
            g1= (s1-k1)/40            
            GFc1 = GF_field[int(round(cell.xCOM)), int(round(cell.yCOM)), int(round(cell.zCOM))]
            cell.targetVolume += varii_1*((g1/30)+(GFc1/20))
            cell.lambdaVolume = 20.0

            if int(round(cell.xCOM))>97 or int(round(cell.xCOM))<3 or int(round(cell.yCOM))>97 or int(round(cell.yCOM))<3:
                self.delete_cell(cell)
                # c+=1
# #                    print('c:',c)
                # if c==1:
                    # mcsOut = mcs
 

        for cell in self.cell_list_by_type(self.CANCER_R):
            neighbor_list_2 = self.get_cell_neighbor_data_list(cell)
            k2 = neighbor_list_2.common_surface_area_with_cell_types(cell_type_list=[1,2])
            s2 = cell.surface
            g2= (s2-k2)/40            
            GFc2 = GF_field[int(round(cell.xCOM)), int(round(cell.yCOM)), int(round(cell.zCOM))]
            cell.targetVolume += varii_2*((g2/30)+(GFc2/20))
            cell.lambdaVolume = 20.0 

            if int(round(cell.xCOM))>97 or int(round(cell.xCOM))<3 or int(round(cell.yCOM))>97 or int(round(cell.yCOM))<3:
                self.delete_cell(cell)
                # c+=1
# #                    print('c:',c)
                # if c==1:
                    # mcsOut = mcs
 
    # ck=c
    # if ck!=0 and mcs%5==0:
        # c_value.append(ck)
# #        print(c_value)
            

    # def finish(self):
        # global mcsOut
        # global c_value
        # self.f=open(CompuCellSetup.getScreenshotDirectoryName() + "\cvalue.txt","w+")
        # #f= open("cvalue.txt","w+")
        # self.f.write("%d\r\n\n"%mcsOut)
        # for ck in c_value:
            # self.f.write("%d\r"%ck)
        # self.f.close()          

#jsp code
    # def finish(self):
        # global mcsOut
        # global c_value
        # # saving number of cells getting deleted at lattice boundary at each mcs
# ##jps## self.f=open(CompuCellSetup.getScreenshotDirectoryName() + "\cvalue.txt","w+")
# ##jps## file_obj, file_path = self.open_file_in_simulation_output_folder("\cvalue.txt", mode='w+')
        # from pathlib import Path
        # if self.output_dir is not None:
            # output_path = Path(self.output_dir).joinpath("cvalue.txt")
            # # create folder to store data
            # output_path.parent.mkdir(parents=True, exist_ok=True)
            # try:
                # file_handle = open(output_path, 'w+')
            # except IOError:
                # print ("Could not open file for writing: ",output_path)
                # return
        
        # print("\n\n\t file_handle: ",file_handle)
        # file_handle.write("%d\r\n\n"%mcsOut) #saving the number of MCS when cells start getting deleted
        # for ck in c_value:
            # file_handle.write("%d\r"%ck)
        # file_handle.close()  




class MatrixDegradation(SteppableBasePy):
    
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
    def start(self):
        pass
        # self.cellA = self.potts.createCell()
        # self.cellA.type = self.A
        # self.cell_field[10:12, 10:12, 0] = self.cellA

        # self.cellB = self.potts.createCell()
        # self.cellB.type = self.B
        # self.cell_field[92:94, 10:12, 0] = self.cellB        
        
        
    def step(self,mcs):
        
        clcell = self.potts.createCell()
        clcell.type = self.C_LYSED
        
        llcell = self.potts.createCell()
        llcell.type = self.L_LYSED        
        
        MMP_field = self.field.MMP
        I_field = self.field.I
        
        global var
        lysed_id = []
        
        lysed_id_2 = []
        
        for cell in self.cell_list_by_type(self.COI):
            cellDict = CompuCell.getPyAttrib(cell)
            MMPc = MMP_field[cell.xCOM, cell.yCOM, cell.zCOM]
            Ic = I_field[cell.xCOM, cell.yCOM, cell.zCOM]
            
            if Ic > 0.0005:
                T1 = (MMPc/Ic)
                
                if T1 > var:
                    cell.type = self.C_LYSED
                    lysed_id.append(cell)
                # #The code below is alternate way to implement matrix degradation with MMP chem field alone with out TIMP.    
                # # else:
            # if MMPc>2.0:
                # cell.type = self.C_LYSED
                # lysed_id.append(cell)
                    
        for cell in lysed_id:
            cellDict = CompuCell.getPyAttrib(cell)
            if hasattr(cell,"mcsL")==True:
                cell.targetVolume -= 0.005
                cell.lambdaVolume = 20.0
                print('Lets see if this part of code is being executed')
            else:
                cellDict = CompuCell.getPyAttrib(cell)
                mcs = self.simulator.getStep()
                mcsL = mcs
                cellDict["mcsL"] = mcsL
                
                
                
        mcs = self.simulator.getStep()
        for cell in self.cell_list:
            cellDict = CompuCell.getPyAttrib(cell)
            if cell.type == self.C_LYSED:
                for val in cellDict.items():
                    cd1 = val[1]
                    
                    if cell.type==self.C_LYSED and mcs == (cd1+40):
                        cell.type = self.NCI
                        


        for cell in self.cell_list_by_type(self.LAMININ):
            cellDict = CompuCell.getPyAttrib(cell)
            MMPc = MMP_field[cell.xCOM, cell.yCOM, cell.zCOM]
            Ic = I_field[cell.xCOM, cell.yCOM, cell.zCOM]
            
            if Ic > 0.0005:
                T1 = (MMPc/Ic)
                
                if T1 > var:
                    cell.type = self.L_LYSED
                    lysed_id_2.append(cell)

        for cell in lysed_id_2:
            cellDict2 = CompuCell.getPyAttrib(cell)
            if hasattr(cell,"mcsL")==True:
                cell.targetVolume -= 0.005
                cell.lambdaVolume = 20.0
                print('Lets see if this part of code is being executed')
            else:
                cellDict2 = CompuCell.getPyAttrib(cell)
                mcs = self.simulator.getStep()
                mcsL = mcs
                cellDict2["mcsL"] = mcsL

                
        mcs = self.simulator.getStep()
        for cell in self.cell_list:
            cellDict2 = CompuCell.getPyAttrib(cell)
            if cell.type == self.L_LYSED:
                for val2 in cellDict2.items():
                    cd2 = val2[1]
                    
                    if cell.type==self.L_LYSED and mcs == (cd2+40):
                        cell.type = self.NCI

                        

# check this below mitosis code with the MITOSIS CC3D DEMOS examples of 4.7.0       
# class MitosisSteppable(MitosisSteppableBase):
    # def __init__(self,frequency=1):
        # MitosisSteppableBase.__init__(self,frequency)

    # def step(self, mcs):

        # cells_to_divide=[]
        # #for cell in self.cell_list:
        # for cell in self.cell_list_by_type(self.CANCER_S, self.CANCER_R):        
            # if cell.volume>30:
                # cells_to_divide.append(cell)

        # for cell in cells_to_divide:

            # self.divide_cell_random_orientation(cell)
            # # Other valid options
            # # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # # self.divide_cell_along_major_axis(cell)
            # # self.divide_cell_along_minor_axis(cell)

    # def update_attributes(self):
        # # reducing parent target volume
        # self.parent_cell.targetVolume /= 2.0                  

        # self.clone_parent_2_child()            

        # # # for more control of what gets copied from parent to child use cloneAttributes function
        # # # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        # # if self.parent_cell.type==1:
            # # self.child_cell.type=2
        # # else:
            # # self.child_cell.type=1



#Two explicit Mitosis steppables
class MitosisSteppable_S(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_S_to_divide=[]
        #for cell in self.cell_list:
        for cell in self.cell_list_by_type(self.CANCER_S):        
            if cell.volume>30:
                cells_S_to_divide.append(cell)

        for cell in cells_S_to_divide:
            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0

        self.clone_parent_2_child()
        # Both parent and child stay CANCER_S: S -> S + S. Mitosis never
        # changes phenotype identity on its own. If a forced-S+S division
        # was pending (from a recent S->R switch elsewhere), this division
        # already satisfies it -- just consume the flag so it doesn't carry
        # over to a later, unrelated division.
        _consume_force_next_division_ss()

        # A newborn cell starts with zero accumulated drug dose.
        # clone_parent_2_child() copied the parent dict, so we must
        # explicitly zero the child's drugUptake (and clear any fate
        # flags the parent may already carry so the child is treated
        # as a fresh, drug-naive cell).
        child_dict = CompuCell.getPyAttrib(self.child_cell)
        child_dict["drugUptake"]  = 0.0
        child_dict.pop("drugHit",     None)
        child_dict.pop("mcsDrugHit",  None)
            
  
class MitosisSteppable_R(MitosisSteppableBase):
    """
    CANCER_R mitosis.

    R -> R + R by default: both parent and child stay CANCER_R after
    division. The only exception: if FORCE_NEXT_DIVISION_SS is pending
    (raised by RandomTypeSwitchSteppable the moment it converts some cell
    S -> R), then THIS division -- whichever R cell happens to divide next,
    anywhere in the lattice -- is forced to yield S + S instead, and the
    parent cell itself is converted to CANCER_S too. The flag is consumed
    (cleared) the instant it's used, so it only overrides a single division.
    """

    # Volume threshold for CANCER_R division (higher → slower division rate)
    CANCER_R_DIVIDE_VOLUME =30

    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_R_to_divide = []
        for cell in self.cell_list_by_type(self.CANCER_R):
            if cell.volume > self.CANCER_R_DIVIDE_VOLUME:
                cells_R_to_divide.append(cell)

        for cell in cells_R_to_divide:
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # Halve parent's target volume after split
        self.parent_cell.targetVolume /= 2.0

        # Clone all attributes (volume params, dict entries) from parent → child
        self.clone_parent_2_child()
        # Parent type stays CANCER_R, and the child is left untouched by
        # clone_parent_2_child() too, so it also stays CANCER_R: R -> R + R
        # -- UNLESS a forced-S+S division is pending (see class docstring).

        if _consume_force_next_division_ss():
            # One-shot override: this division must yield S + S, even
            # though the dividing cell was CANCER_R. Convert both parent
            # and child to CANCER_S.
            self.parent_cell.type = self.CANCER_S
            self.child_cell.type  = self.CANCER_S
            self.parent_cell.targetSurface = 20.0
            self.parent_cell.lambdaSurface = 2.0
            self.child_cell.targetSurface  = 20.0
            self.child_cell.lambdaSurface  = 2.0

        # Reset the child's uptake counter so it begins accumulating drug
        # independently from its parent's history.
        child_dict = CompuCell.getPyAttrib(self.child_cell)
        child_dict["drugUptake"] = 0.0
        child_dict.pop("drugHit",    None)
        child_dict.pop("mcsDrugHit", None)



#Motility steppable
class CellMotilitySteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        # Make sure ExternalPotential plugin is loaded
        # negative lambdaVecX makes force point in the positive direction

        for cell in self.cell_list_by_type(self.CANCER_S, self.CANCER_R):        
            # force component pointing along X axis
            cell.lambdaVecX = 10.1 * uniform(-1.0, 1.0)
            # force component pointing along Y axis
            cell.lambdaVecY = 10.1 * uniform(-1.0, 1.0)


        
######################### Secretion steppable ##########################
class SecretionSteppable(SecretionBasePy):
    
    def __init__(self,frequency=1):
        SecretionBasePy.__init__(self,frequency)
        
    def start(self):
        pass
        # self.fieldNameMMP='MMP'
        # self.fieldNameI='I'

    def step(self,mcs):
        
        MMPsecretor  = self.get_field_secretor("MMP")
        Isecretor    = self.get_field_secretor("I")
        GFsecretor   = self.get_field_secretor("GF")
        DRUGsecretor = self.get_field_secretor("DRUG")
        
        MMP_field  = self.field.MMP
        I_field    = self.field.I
        DRUG_field = self.field.DRUG

        for cell in self.cell_list:

            MMPc = MMP_field[cell.xCOM, cell.yCOM, cell.zCOM]
            Ic   = I_field[cell.xCOM, cell.yCOM, cell.zCOM]
            
            if cell.type == self.CANCER_S:
                
                x = random.randint(0,4)
                A1 = 0.25
                I1 = A1

                MMPsecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,A1,[self.COI])
                Isecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,I1,[self.COI])
                MMPsecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,A1,[self.LAMININ])
                Isecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,I1,[self.LAMININ])
                GFsecretor.uptakeInsideCell(cell,0.1,0.1)

                # ── Drug uptake for CANCER_S ──────────────────────────────────
                # Read local drug concentration BEFORE uptake, then internalise
                # a fraction of it (this physically removes drug from the field).
                # The absorbed amount is accumulated in cell_dict["drugUptake"],
                # which DrugResponseSteppable then reads to trigger fate decisions
                # instead of the raw instantaneous field value.
                drug_local_S = DRUG_field[int(round(cell.xCOM)),
                                          int(round(cell.yCOM)),
                                          int(round(cell.zCOM))]
                DRUGsecretor.uptakeInsideCell(cell, DRUG_UPTAKE_RATE_S, 1.0)
                cell_dict_S = CompuCell.getPyAttrib(cell)
                cell_dict_S["drugUptake"] = (cell_dict_S.get("drugUptake", 0.0)
                                             + DRUG_UPTAKE_RATE_S * drug_local_S)

            if cell.type == self.CANCER_R:

                x = random.randint(0,4)
                A2 = 0.25
                I2 = A2

                MMPsecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,A2,[self.COI])
                Isecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,I2,[self.COI])
                MMPsecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,A2,[self.LAMININ])
                Isecretor.secreteOutsideCellAtBoundaryOnContactWith(cell,I2,[self.LAMININ])
                GFsecretor.uptakeInsideCell(cell,0.1,0.1)

                # ── Drug uptake for CANCER_R (5× lower — resistance phenotype) ─
                # CANCER_R internalises much less drug per MCS.  The accumulator
                # is still tracked in cell_dict["drugUptake"] for observability /
                # future use, but DrugResponseSteppable does not act on R cells.
                drug_local_R = DRUG_field[int(round(cell.xCOM)),
                                          int(round(cell.yCOM)),
                                          int(round(cell.zCOM))]
                DRUGsecretor.uptakeInsideCell(cell, DRUG_UPTAKE_RATE_R, 1.0)
                cell_dict_R = CompuCell.getPyAttrib(cell)
                cell_dict_R["drugUptake"] = (cell_dict_R.get("drugUptake", 0.0)
                                             + DRUG_UPTAKE_RATE_R * drug_local_R)   
   
 
            if cell.type==self.L_LYSED:
                GFsecretor.secreteInsideCell(cell,0.5)
                MMPsecretor.uptakeInsideCell(cell,1.0,1.0)
                Isecretor.uptakeInsideCell(cell,0.5,1.0)  

            if cell.type==self.C_LYSED:                
                GFsecretor.secreteInsideCellAtBoundary(cell,1.0)
                MMPsecretor.uptakeInsideCell(cell,1.5,1.0)
                Isecretor.uptakeInsideCell(cell,0.5,1.0)
                




class OrientedConstraintSteppable(SteppableBasePy):
    def __init__(self, frequency):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list_by_type(self.COI):
            o = random.randint(2,10)
            print(o)
            cell.lambdaVolume = 20.0
            cell.targetVolume = cell.volume

            # Here, we define the axis of elongatino.
            self.orientedGrowthPlugin.setElongationAxis(cell, math.cos(math.pi/o), math.sin(math.pi/o))
            # And this function gives a 2 pixel width to each cell
            self.orientedGrowthPlugin.setConstraintWidth(cell, 1.0)

            # Make sure to enable or disable elongation in all cells
            self.orientedGrowthPlugin.setElongationEnabled(cell, True)



######################### Drug response steppable ##########################
class DrugResponseSteppable(SteppableBasePy):
    """
    Models the effect of accumulated drug uptake on CANCER_S cells.

    Instead of reacting to the instantaneous local DRUG field concentration,
    each CANCER_S cell now carries a per-cell accumulator: cell_dict["drugUptake"].
    SecretionSteppable increments this every MCS by (DRUG_UPTAKE_RATE_S × local [DRUG]).

    Fate decision (triggered once, when drugUptake crosses DRUG_UPTAKE_THRESHOLD):
      • the cell is deterministically condemned ('drugHit') and loses
        targetVolume by DRUG_DEATH_RATE each MCS until it is deleted.
      • there is no longer a stochastic branch to CANCER_R here -- S/R
        identity is now governed *only* by RandomTypeSwitchSteppable and by
        mitosis (S -> S+S, R -> R+R), never by the drug response itself.

    CANCER_R cells are intentionally left untouched (they do accumulate
    drugUptake for observability, but no fate decision is applied).
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Initialise drugUptake=0 on all cancer cells that already exist
        for cell in self.cell_list_by_type(self.CANCER_S, self.CANCER_R):
            cell_dict = CompuCell.getPyAttrib(cell)
            cell_dict.setdefault("drugUptake", 0.0)

    def step(self, mcs):

        cells_to_delete = []

        for cell in self.cell_list_by_type(self.CANCER_S):
            cell_dict = CompuCell.getPyAttrib(cell)

            # Ensure the key exists even if SecretionSteppable hasn't run yet
            cell_dict.setdefault("drugUptake", 0.0)

            # ── One-time fate decision when cumulative uptake crosses threshold ──
            # Deterministic: crossing the threshold always condemns the cell
            # to death. There is no more stochastic resistance branch.
            if not cell_dict.get("drugHit", False):
                if cell_dict["drugUptake"] > DRUG_UPTAKE_THRESHOLD:
                    cell_dict["drugHit"]    = True
                    cell_dict["mcsDrugHit"] = mcs

            # ── Progressive volume loss for cells condemned to die ──
            if cell_dict.get("drugHit", False):
                cell.targetVolume = max(cell.volume - DRUG_DEATH_RATE, 0.0)
                cell.lambdaVolume = 20.0

                if cell.volume < DRUG_MIN_VOLUME:
                    cells_to_delete.append(cell)

        for cell in cells_to_delete:
            self.delete_cell(cell)

        # CANCER_R (resistant) cells are intentionally left untouched here
        

class AreaTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Create the plot window
        self.plot_win = self.add_new_plot_window(
            title='Total Cell Area Over Time',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Total Area (pixels)',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True
        )
        # Add a data series - pass empty lists for x and y data
        self.plot_win.add_plot("TotalArea",color="red",style="lines",size=5)

    def step(self, mcs):
        total_area = 0
        for cell in self.cell_list:
            total_area += cell.volume

        self.plot_win.add_data_point("TotalArea", mcs, total_area)


######################### Random type-switch steppable ##########################
class RandomTypeSwitchSteppable(SteppableBasePy):
    """
    Periodically flips the phenotype (CANCER_S <-> CANCER_R) of a randomly
    chosen cell from the current cancer-cell population.

    This is completely independent of drug exposure and of mitosis -- it is
    the only mechanism in the model that ever changes a cell's S/R identity
    directly:
        • DrugResponseSteppable only kills CANCER_S cells past threshold.
        • Mitosis is normally S -> S+S and R -> R+R, EXCEPT that the one
          division immediately following an S -> R conversion made here is
          forced to produce S+S regardless of which cell divides (see
          FORCE_NEXT_DIVISION_SS / _consume_force_next_division_ss()).

    Timing: the waiting time between successive switch events is drawn from
    Normal(MU_SWITCH_WINDOW, SIGMA_SWITCH_WINDOW) (MCS), and is redrawn after
    every event, so events land at irregular intervals rather than on a
    fixed period, but average out to about once every MU_SWITCH_WINDOW mcs.

    Target: the cell to flip is drawn uniformly at random directly from the
    current CANCER_S + CANCER_R population (no spatial coordinate lookup),
    so every live cancer cell has an equal chance of being picked regardless
    of where it sits on the lattice.
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.next_switch_mcs = None

    def _draw_next_window(self):
        window = random.normalvariate(MU_SWITCH_WINDOW, SIGMA_SWITCH_WINDOW)
        return max(MIN_SWITCH_WINDOW, window)

    def start(self):
        self.next_switch_mcs = self._draw_next_window()

    def step(self, mcs):
        if mcs < self.next_switch_mcs:
            return

        # Pick a random cell directly from the current cancer population.
        candidates = list(self.cell_list_by_type(self.CANCER_S, self.CANCER_R))

        if candidates:
            cell = random.choice(candidates)
            if cell.type == self.CANCER_S:
                cell.type = self.CANCER_R
                # This S -> R conversion means the *next* division anywhere
                # in the lattice -- regardless of which cell's turn it is --
                # must yield S + S instead of the usual rule for that parent.
                global FORCE_NEXT_DIVISION_SS
                FORCE_NEXT_DIVISION_SS = True
            else:
                cell.type = self.CANCER_S

            # Keep surface constraints consistent after the type flip
            # (both types use the same targetSurface/lambdaSurface anyway).
            cell.targetSurface = 20.0
            cell.lambdaSurface = 2.0

        # Schedule the next event regardless of whether any cancer cells
        # currently exist (candidates could be empty, e.g. early in the run
        # or if the population dies out).
        self.next_switch_mcs = mcs + self._draw_next_window()