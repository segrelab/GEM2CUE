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
    | strains <list>: List of Strain objects
    | prep_kwargs <dict>: Keyword arguments to pass to `prep_cobrapy_models()` from surfinFBA
        * See bottom of: https://github.com/jdbrunner/surfin_fba
    """
    def __init__(self, strains: List[Strain], init_strain_biomasses: Union[list, dict], ):
        # Check that an 'e' compartment is in the model
        for s in strains:
            if 'e' not in s.model.compartments:
                raise KeyError('An extracellular compartment, "e" must be defined in each model.')
        
        self.strains = strains
        self.strain_dict = {s.name: s.model for s in self.strains}