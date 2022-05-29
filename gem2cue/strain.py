import cobra
import numpy as np

class Strain:
    """
    Object for single strain

    | name <str>: Name of strain
    | model <cobra.Model>: Cobrapy model associated with strain. The extracellular compartment must be labeled 'e' for dFBA simulations.
    | metabolite_uptake <dict>: Uptake rate of each exchange reaction in `model.medium` for dFBA simulations.
        * If provided, uptake reacts for all reactions in `model.medium` must be specified
        
        Ex. {'EX_nh4_e': 2.3,
            'EX_co2_e': 1,
            'EX_glc__D_e': 1,
            'EX_h_e': 3,
            'EX_h2o_e': 2,
            'EX_o2_e': .5,
            'EX_pi_e': 1}
        
        * If not provided, defaults to random values between 0 and 1
    """
    def __init__(self, name: str, model: cobra.Model, metabolite_uptake: dict=None, gc_content: float=None, genome_length: int=None):
        self.name = name
        self.model = model
        self.gc_content = gc_content
        self.genome_length = genome_length

        # Handle metabolite uptake
        # If not defined, assign random numbers 
        #   (like surfinFBA does: https://github.com/jdbrunner/surfin_fba/blob/1566282ddb628be3914e54b6ccd4468958338699/surfinFBA/Surfin_FBA.py#L350)
        #   TODO: Metabolite uptakes should be same across strains (like surfinFBA does, except they also only take the intersection across strains...)
        if not metabolite_uptake:
            metabolite_uptake = np.random.rand(len(self.model.medium))
        
        # If missing reactions, through error
        if metabolite_uptake.keys() != self.model.medium.keys():
            raise KeyError('If providing uptake rates, uptake for all metabolites in `model.medium` must be specified')