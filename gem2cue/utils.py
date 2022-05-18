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


class Organism:
    "A metabolic model for a single organism in a Community"

    def __init__(self, strain: Strain, biomass: float = 0.1, flux_parameters: dict = None):
        """
        strain: A Strain class object with a cobrapy model 

        biomass: Initial biomass
        flux_parameters: Dictionary of vmax and km (default: {'vmax': 2.0, 'km': 0.5})
        """
        self.strain = strain
        self._exchange_reaction_ids = [e.id for e in model.exchanges]
        self._biomasses = [biomass]
        if not flux_parameters:
            flux_parameters = {'vmax': 2.0, 'km': 0.5}
        self.flux_parameters = flux_parameters

        self.growth_rates = []
        self._fluxes = []
        self._uptake_concentrations = []

        self.infeasible_timestep = None

    @property
    def biomass(self) -> float:
        return self._biomasses[-1]
    
    @property
    def biomasses(self) -> pd.Series:
        return pd.Series(self._biomasses).rename(self.model.id)
    
    @property
    def fluxes(self) -> pd.DataFrame:
        fluxes_df = pd.concat(self._fluxes, axis=1).T.reset_index(drop=True).rename_axis(index='timestep', columns='reaction')
        non_zero_fluxes = fluxes_df.ne(0).any().pipe(lambda x: x[x]).index
        # Subset to non-zero fluxes
        return fluxes_df[non_zero_fluxes]

    @property
    def uptake_concentrations(self) -> pd.DataFrame:
        return -pd.DataFrame(self._uptake_concentrations).rename_axis(index='timestep', columns='reaction')

class Experiment:
    "A collection of metabolic model(s) in an environment"

    def __init__(self, organisms: List[Organism], media: Media, timepoints: int = 20, dt: float = 0.1):
        """
        organisms: A list of Organism(s)
        media: A `Media` object
        timepoints: Number of timepoints in dFBA simulation
        dt: Size of timestep for each timepoint
        flux_parameters: Dictionary of vmax and km (default: {'vmax': 2.0, 'km': 0.5})
        """
        if isinstance(organisms, Organism):
            organisms = [organisms]
        self.organisms = organisms
        self.Media = media
        self.timepoints = timepoints
        self.dt = dt

        self._timestep = 0
    
    @property
    def biomasses(self):
        "Biomass of each Organism at each timestep"
        return pd.concat([o.biomasses for o in self.organisms], axis=1).rename_axis(index='timestep', columns='organism')

    def _dFBAstep(self, organism: Organism):
        "A single step of dFBA for a single Organism including updates to metabolite concentrations in media"

        "Update metabolite uptake rate given environment concentration"
        # Calculate MM uptake rates for limiting nutrients
        for reaction, concentration in self.Media.media.items():
            if np.isfinite(concentration) & (reaction in organism._exchange_reaction_ids):
                # Compute MM uptake rate for limiting nutrient
                uptake_rate = organism.flux_parameters['vmax'] * concentration / (organism.flux_parameters['km'] + concentration)
                # Impose rate as maximum uptake rate
                organism.model.exchanges.get_by_id(reaction).lower_bound = -uptake_rate

        "Compute SS FBA for this timepoint"
        # Run FBA
        solution = organism.model.optimize()
        if solution.status == 'optimal':
            # Growth rate
            growth = solution.objective_value
            organism.growth_rates.append(growth)
            # Exchange fluxes
            fluxes_now = solution.fluxes[organism._exchange_reaction_ids]
            organism._fluxes.append(fluxes_now)
        else:  
            # If infeasible
            if not organism.infeasible_timestep:
                organism.infeasible_timestep = self._timestep
            growth = 0
            fluxes_now = pd.Series({m: 0 for m in self.Media.media})

        "Update biomass and media concentrations"
        # Update uptake concentrations in organism
        uptake_concentrations_now = {}
        for m, flux in fluxes_now.items():
            if flux != 0:
                # Calculate uptake concentration in Organism
                uptake_concentration = flux * organism.biomass * self.dt
                uptake_concentrations_now[m] = uptake_concentration
        organism._uptake_concentrations.append(uptake_concentrations_now)
        # Uptake concentrations in media
        self.Media.useMedia(uptake_concentrations_now)

        # Update biomass
        biomass = organism.biomass + growth * organism.biomass * self.dt
        organism._biomasses.append(biomass)

    def dFBA(self):
        "Run dFBA for all model(s)"
        for t in range(self.timepoints):
            self._timestep = t
            for organism in self.organisms:
                self._dFBAstep(organism)