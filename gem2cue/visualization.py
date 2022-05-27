import matplotlib as plt
from matplotlib.sankey import Sankey

def sankey(self, co2_rxn: str ='EX_co2_e', biomass_rxn: str = 'BIOMASS_Ecoli_core_w_GAM'):
        """ Create and save a sankey diagram for the CUE of a specific model
    
        Args:
        model (cobra.core.Model): Model to be solved
        co2_rxn (str): Name of the respiration reaction in the model (defaults to 
            'EX_co2_e')
        biomass_rxn (str): Name of the biomass reaction in the model (defaults to
            'BIOMASS_Ecoli_core_w_GAM')
        
        Returns:
        plt object
        """
        # If the experiment does not have a solution, run it
        if self.solution is None:
            self.run()

        # Look up the exchange reactions
        ex_c_atoms = self._atomExchangeMetabolite()

        # Find the uptake flux
        uptake_rxns = self.solution.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x<0]).index # Negative flux for uptake
        uptake = abs(sum([self.solution.get_primal_by_id(r) * ex_c_atoms[r] for r in uptake_rxns if r != co2_rxn])) # Doesn't really matter if we say no CO2 rxn or not, it shouldn't be a negative flux

        # Find the respiration flux
        resp = -1 * self.solution.get_primal_by_id(co2_rxn) # Don't multiply because we know that the number of c atoms is 1 for CO2

        # Find the exudation flux
        exudation_rxns = self.solution.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x>0]).index # Postive flux for exudation
        ex = -1 * abs(sum([self.solution.get_primal_by_id(r) * ex_c_atoms[r] for r in exudation_rxns if r != co2_rxn])) # Do really need to say no CO2 reaction

        # Find the biomass flux
        # Make a dictionary of the fluxes for each metabolite in the biomass reaction
        biomass_dict = getattr(self.strain.model.reactions, biomass_rxn).metabolites
        # Find the number of C atoms that the biomass reaction takes in
        biom_in = abs(sum([biomass_dict[m] * self.solution.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] < 0 if 'C' in m.elements]))
        # Find the number of C atoms that the biomass reaction puts out
        biom_out = abs(sum([biomass_dict[m] * self.solution.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] > 0 if 'C' in m.elements]))
        # Find the difference
        # Is this what actually goes to biomass?
        biom = -1 * (biom_in - biom_out)

        # Normalize to uptake
        uptake_norm = 1
        resp_norm = resp / uptake
        ex_norm = ex / uptake
        biom_norm = biom / uptake

        # Make the plot
        fig = plt.figure()
        Sankey(flows = [uptake_norm, resp_norm, ex_norm, biom_norm],
            labels = ['U = A', 'R', 'EX', 'G'],
            orientations = [0, 1, 1, 0],
            pathlengths = [1, 0.25, 0.25, 0.5]).finish()
        plt.title('Carbon Flow for ' + self.strain.name)
        return(fig)