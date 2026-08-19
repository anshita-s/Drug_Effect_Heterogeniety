import sys
from os import environ
from cc3d import CompuCellSetup
from cc3d.cpp import CompuCell
import string
from os import getcwd
from drugSteppables import CellLayoutSteppable
CompuCellSetup.register_steppable(steppable=CellLayoutSteppable(frequency=1))
from drugSteppables import VolumeParamSteppable
CompuCellSetup.register_steppable(steppable=VolumeParamSteppable(frequency=1))
from drugSteppables import DrugResponseSteppable
CompuCellSetup.register_steppable(steppable=DrugResponseSteppable(frequency=1))
from drugSteppables import MatrixDegradation
CompuCellSetup.register_steppable(steppable=MatrixDegradation(frequency=1))
from drugSteppables import MitosisSteppable_S
CompuCellSetup.register_steppable(steppable=MitosisSteppable_S(frequency=1))
from drugSteppables import MitosisSteppable_R
CompuCellSetup.register_steppable(steppable=MitosisSteppable_R(frequency=1))
from drugSteppables import CellMotilitySteppable
CompuCellSetup.register_steppable(steppable=CellMotilitySteppable(frequency=10))
from drugSteppables import SecretionSteppable
CompuCellSetup.register_steppable(steppable=SecretionSteppable(frequency=4))
from drugSteppables import OrientedConstraintSteppable
CompuCellSetup.register_steppable(steppable=OrientedConstraintSteppable(frequency=1))
from drugSteppables import AreaTrackerSteppable
CompuCellSetup.register_steppable(steppable=AreaTrackerSteppable(frequency=1))
from drugSteppables import RandomTypeSwitchSteppable
CompuCellSetup.register_steppable(steppable=RandomTypeSwitchSteppable(frequency=1))
CompuCellSetup.run()