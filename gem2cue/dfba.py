from ctypes import Union
from strain import Strain
from typing import List


try:
    import surfinFBA
except:
    msg = """surfinFBA is required for dFBA simulations.
    To install surfinFBA:
        1) Clone the repository: git clone git@github.com:jdbrunner/surfin_fba.git
        2) Enter the directory: cd surfin_fba
        3) Install: pip install .

    surfinFBA requires either CPLEX or GUROBI. 
    """
    raise ImportError(msg)


class Timecourse:
    """
    Modeling CUE of a strain or community over time with surfinFBA

    Inputs:
    | strains <list>: List of Strains
    | x0 <list or dict>: Initial biomass for each Strain
    | y0
    """
    def __init__(self, strains: List[Strain], x0: dict, endtime, y0: dict=None, met_in=None, met_out=None, metabolite_names=[], report_activity=False, detail_activity=False, initres=0.001, concurrent=True, solver='both',  enoughalready=10, flobj=None, chk_round=5):
        self.model_dict = {s.name: s.model for s in strains}
