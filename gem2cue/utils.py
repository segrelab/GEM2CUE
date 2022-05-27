"Class objects for running dFBA- copied from Michael's dFBA package"

from typing import List
import cobra
import numpy as np
import pandas as pd
import warnings
import matplotlib as plt
from matplotlib.sankey import Sankey

class Media:
    "Environmental media for a Community over time"

    def __init__(self, media: dict = None):
        """
        media: A dictionary of cobrapy exchange reactions
        """
        self.media = media


class Strain:
    "A model and it's associated metadata"

    def __init__(self, name: str, model: cobra.core.Model, metadata: dict = None):
        """
        name:
        model:
        gc_content:
        genome_length:
        """
        self.name = name
        self.model = model.copy()
        self.metadata = metadata

    def update_medium(self, new_medium: Media):
        clean_media = {i: v for i, v in new_medium.media.items() if i in self.model.exchanges}
        self.model.medium = clean_media


class Experiment:
    "A collection of one strain in an environment"

    def __init__(self, strain: Strain):
        """
        organisms: A list of Organism(s)
        solution
        cue
        """
        self.strain = strain
        self.solution = None
        self.cue = None

    def run(self):
        "Run FBA"
        # Warn if the experiment already has a solution
        if self.solution is not None:
            warnings.warn('There is already a solution saved to this experiment, running will overwrite those results')

        # Solve FBA
        sol = self.strain.model.optimize()

        # Update the experiment object
        self.solution = sol

    def _atomExchangeMetabolite(self, atom = 'C', ex_nomenclature = {'e'}):
        # TODO: Infer the ex_nomenclature rather than forcing the user to provide it
        """Get number of carbon atoms associated with each exchange reaction
    
        Args:
            model (cobra.core.Model): A file that has already been read in
            atom (string): String of atom of interest

        Returns:
            ex_atoms (dict): Dictionary with the IDs as the rxn names, and the
                values as the number of atoms associates
        """
        # FIXME: This is where the issue is
        # compartment for CarveMe models is C_e
        # Compartment for BiGG models is e
        ex_atoms = {r.id: m.elements[atom] for m in self.strain.model.metabolites for r in m.reactions if atom in m.elements if r.compartments == ex_nomenclature}
        
        return ex_atoms

    def CUE(self, co2_rxn='EX_co2_e', ex_nomenclature = {'e'}, definition = 'rCUE'):
        """ Calculate CUE, rCUE uses assumes that respiration is the only waste,
        using the formula: CUE = 1 - CO2 / Uptake C. Any other definition will
        compute CUE as Gross Growth Efficiency using the following definition

                    sum(uptake C) - sum(secretion C)
            CUE =  ----------------------------------
                            sum(uptake C)

        Args:
        self (Experiment): Experiment with the model/media to use
        co2_rxn (string): Name of the respiration reaction in the model
        ex_nomenclature (string): Nomenclature the model uses to denote exchange
            reactions (BiGG uses 'e', CarveMe uses 'C_e')
        definition (string): What definition to use
        """
        # Warn if the experiment already has a CUE value
        if self.cue is not None:
            warnings.warn('There is already a cue value saved to this experiment, running will overwrite those results')

        # If the experiment does not have a solution, run it
        if self.solution is None:
            self.run()

        # Get C atoms for each exchange reaction
        ex_c_atoms = self._atomExchangeMetabolite(ex_nomenclature = ex_nomenclature)

        # Calculate CUE using one of the following definitions
        if definition == 'rCUE':
            # Subset to uptake reactions (negative flux)
            uptake_rxns = self.solution.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x<0]).index
            # Calculate uptake C flux
            uptake = sum([self.solution.get_primal_by_id(r) * ex_c_atoms[r] for r in uptake_rxns if r != co2_rxn])
            if uptake == 0:
                cue = None
            else:
                # Calculate CUE
                cue = 1 - abs(self.solution.get_primal_by_id(co2_rxn) / uptake)
        else:
            # Assume that the only other option is GGE
            # Get C fluxes (flip signs so that uptake is positive)
            c_ex_fluxes = np.array([self.solution.get_primal_by_id(r) * -c for r, c in ex_c_atoms.items()])
            # Calculate GGE
            cue = c_ex_fluxes.sum() / c_ex_fluxes[c_ex_fluxes > 0].sum()

        # Update the experiment with the results
        self.cue = cue


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