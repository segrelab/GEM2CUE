"Class objects for running dFBA- copied from Michael's dFBA package"

from typing import List
import cobra
import numpy as np
import pandas as pd


class Media:
    "Environmental media for a Community over time"

    def __init__(self, media: dict = None, fixed: List[str] = None):
        """
        media: A dictionary of cobrapy exchange reactions and starting concentration. 
        | * LIMITING NUTRIENTS: Provide desired concentration
        | * NON-LIMITING NUTRIENTS: Indicate with unlimited supply with `np.inf` and completely limited with `-np.inf`
         Ex: Limited glucose with unlimited CO2, H+, H2O, NH4, O2, and Pi:
            media = {'EX_glc__D_e': 10.0, 'EX_co2_e': inf, 'EX_h_e': inf, 'EX_h2o_e': inf, 'EX_nh4_e': inf, 'EX_o2_e': inf, 'EX_pi_e': inf}
        fixed: A list of metabolites with fixed concentration (ex. low O2)
        """

        self._medias = [media]
        if not fixed:
            fixed = []
        self._fixed = fixed

    @property
    def media(self):
        return self._medias[-1]

    @property
    def medias(self):
        return pd.DataFrame(self._medias).rename_axis(index='timestep', columns='reaction')

    def useMedia(self, uptakes: dict):
        "Adjust media by `uptakes`, a dictionary of metabolites and uptake concentrations"
        new_concentrations = self.media.copy()

        for m, u in uptakes.items():
            if m in self._fixed:
                new_c = self.media[m]
            else:
                new_c = self.media.get(m, 0) + u
                if new_c < 0:
                    new_c = 0
            new_concentrations[m] = new_c
        
        # Append to `medias` list
        self._medias.append(new_concentrations)

class Strain:
    "A model and it's associated metadata"

    def __init__(self, name, model, gc_content, genome_length):
        """
        name:
        model:
        gc_content:
        genome_length:
        """
        self.name = name
        self.model = model
        self.gc_content = gc_content
        self.genome_length = genome_length


class Experiment:
    "A collection of one strain in an environment"

    def __init__(self, strain: Strain, media: Media):
        """
        organisms: A list of Organism(s)
        media: A `Media` object
        """
        self.strain = strain
        self.Media = media


    def atomExchangeMetabolite(self, atom = 'C'):
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
        ex_atoms = {r.id: m.elements[atom] for m in self.strain.model.metabolites for r in m.reactions if atom in m.elements if r.compartments == {'e'}}
        
        return ex_atoms