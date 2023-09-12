import pandas as pd
import cobra

def get_c_ex_rxns(model, atom = 'C'):
    # TODO: Infer the ex_nomenclature rather than forcing the user to provide it
    """Get number of carbon atoms associated with each exchange reaction

    Args:
        model (cobra.core.Model): A file that has already been read in
        atom (string): String of atom of interest
        ex_nomenclature (set): Set of strings of the exchange nomenclature

    Returns:
        ex_atoms (dict): Dictionary with the IDs as the rxn names, and the
            values as the number of atoms associates
    """
    # This assumes that there is only ever one reactant in an exchange
    # reaction
    ex_atoms = {ex_rxn.id: ex_rxn.reactants[0].elements[atom] for 
                ex_rxn in model.exchanges
                if atom in ex_rxn.reactants[0].elements}

    return ex_atoms


def get_c_ex_rxn_fluxes(solution, c_ex_rxns, tool_used):
    """Use the list of dictionary of exchange reactions to make a
    two dictionaries of reactions and the number of carbon atoms
    exchanged (the flux times the number of carbon atoms in the
    metabolite). Makes one dictionary for uptake exchanges and another
    dictionary for secretion excahnges.

    Args:
    solution (pd.Series OR cobra.Solution): Results from FBA
    tool_used (str): Which tool was used to run FBA. Options are
        "COMETS" or "COBRApy" (Capitalization does not matter).

    Returns:
    uptake (dict):
    secretion (dict):
    
    """
    # Check that the tool used matches with the expected type of the
    # solution object
    if tool_used.lower() != 'comets' and tool_used.lower() != 'cobrapy':
        raise ValueError('Function does not recognize the value supplied ' +
                         'for `tool_used`. Select from `COMETS` or `COBRApy`' +
                         ' (capitalization does not matter). You supplied ' +
                         tool_used + '.')
    if tool_used.lower() == 'comets' and not isinstance(solution, pd.Series):
        raise ValueError('Function was expecting results from a COMETS' +
                         'simulation as a pandas.Series object, but was ' +
                         'given the solution as a ' + type(solution) +
                         'object.')
    if tool_used.lower() == 'cobrapy' and not isinstance(solution, cobra.Solution):
        raise ValueError('Function was expecting results from a COBRApy' +
                         'simulation as a cobra.Solution object, but was ' +
                         'given the solution as a ' + type(solution) +
                         'object.')
    
    # Make a dictionary of all carbon atom fluxes (reaction flux * number of
    # carbon atoms). Have to access the fluxes differently if the results are
    # from COMETS or COBRApy.
    if tool_used.lower() == 'comets':
        c_ex_fluxes = {rxn_id: float(solution[rxn_id]) * n_atoms
                       for rxn_id, n_atoms in c_ex_rxns.items()}
    if tool_used.lower() == 'cobrapy':
        c_ex_fluxes = {rxn_id: solution.fluxes[rxn_id] * n_atoms
                       for rxn_id, n_atoms in c_ex_rxns.items()}
    
    # Set the signs for the uptake and secretion
    
    # Separate the uptake and secretion fluxes based on the sign
    uptake = {rxn_id: abs(atom_flux)
              for rxn_id, atom_flux in c_ex_fluxes.items()
              if atom_flux < 0}
    secretion = {rxn_id: abs(atom_flux)
                 for rxn_id, atom_flux in c_ex_fluxes.items()
                 if atom_flux > 0}

    return uptake, secretion

def get_co2_secretion(secretion_fluxes, co2_ex_rxn = 'EX_co2_e'):
    """Get the total number of carbon atoms lost from the cell as CO2

    Args:
    secretion_fluxes (dict): Dictionary of carbon secreting reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    co2_ex_rxn (str): Reaction ID for the CO2 exchange reaction.
        Defaults to the BiGG ID 'EX_co2_e'.

    Returns:
    co2_flux (float): Numeric value for the carbon atom flux for CO2
        secretion
    
    """
    if co2_ex_rxn not in secretion_fluxes.keys():
        co2_flux = 0
    else:
        co2_flux = secretion_fluxes[co2_ex_rxn]

    return co2_flux


def get_org_c_secretion(secretion_fluxes, co2_ex_rxn = 'EX_co2_e'):
    """Get the total number of carbon atoms lost from the cell as
    organic carbon

    Args:
    secretion_fluxes (dict): Dictionary of carbon secreting reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    co2_ex_rxn (str): Reaction ID for the CO2 exchange reaction.
        Defaults to the BiGG ID 'EX_co2_e'.

    Returns:
    org_c_secretion_flux (float): Numeric value for the carbon atom flux
        for all secretion reactions other than CO2.
    """
    org_c_secretion_flux = sum([c_atom_flux
                                for rxn, c_atom_flux in secretion_fluxes.items()
                                if rxn != co2_ex_rxn])

    return org_c_secretion_flux


def get_c_uptake(uptake_fluxes):
    """Get the total number of organic carbon atoms taken up. Would
    include CO2 if that is taken up.

    Args:
    uptake_fluxes (dict): Dictionary of carbon importing reactions
        with the reaction ID and the absolute value of the carbon atom
        flux

    Returns:
    uptake_flux (float): Numeric value for the total carbon atom flux
        for all import reactions
    """
    uptake_flux = sum([c_atom_flux
                       for rxn, c_atom_flux in uptake_fluxes.items()])

    return uptake_flux


def get_biomass_carbon(solution, biomass_rxn, model, tool_used = 'COMETS'):
    """Get the total number of carbon atoms used by the biomass reaction

    Args:
    solution (pd.Series OR cobra.Solution): Results from FBA
    biomass_rxn (str): Reaction ID for the biomass reaction
    model (cobra.Model): COBRA model used
    tool_used (str): Which tool was used to run FBA. Options are
        "COMETS" or "COBRApy" (Capitalization does not matter).

    Returns:
    
    """
    # Check that the tool used matches with the expected type of the
    # solution object
    if tool_used.lower() != 'comets' and tool_used.lower() != 'cobrapy':
        raise ValueError('Function does not recognize the value supplied ' +
                         'for `tool_used`. Select from `COMETS` or `COBRApy`' +
                         ' (capitalization does not matter). You supplied ' +
                         tool_used + '.')
    if tool_used.lower() == 'comets' and not isinstance(solution, pd.Series):
        raise ValueError('Function was expecting results from a COMETS' +
                         'simulation as a pandas.Series object, but was ' +
                         'given the solution as a ' + type(solution) +
                         'object.')
    if tool_used.lower() == 'cobrapy' and not isinstance(solution, cobra.Solution):
        raise ValueError('Function was expecting results from a COBRApy' +
                         'simulation as a cobra.Solution object, but was ' +
                         'given the solution as a ' + type(solution) +
                         'object.')

    # Get the flux through the biomass reaction
    if tool_used.lower() == 'comets':
        rxn_flux = solution[biomass_rxn]
    if tool_used.lower() == 'cobrapy':
        rxn_flux = solution.fluxes[biomass_rxn]

    # Get the actual reaction object for the biomass reaction
    rxn_obj = model.reactions.get_by_id(biomass_rxn)

    c_atom_flux = 0
    # Loop through all of the biomass components
    for component, s_coeff in rxn_obj.metabolites.items():
        print('--------------\n' + component.name + '\nCoeff: ' + str(s_coeff))
        # If the component does not contain carbon, skip it
        if 'C' not in component.elements.keys():
            continue
        # Get the number of carbon atoms in the component
        n_c_atoms = component.elements['C']
        print('num. C atoms: ' + str(n_c_atoms))
        # Multiply the number of carbon atoms by the stoichiometric coefficient
        component_flux = n_c_atoms * s_coeff
        # Add the flux to the total c_atom_flux
        c_atom_flux += component_flux
        print('new atom flux: ' + str(c_atom_flux))

    # The final c atom flux is the 
    return abs(c_atom_flux * rxn_flux)


def calculate_cue(uptake_fluxes, secretion_fluxes, co2_ex_rxn = 'EX_co2_e'):
    """Calculate the CUE by using the uptake and secretion dictionaries

    Args:
    uptake_fluxes (dict): Dictionary of carbon importing reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    secretion_fluxes (dict): Dictionary of carbon secreting reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    co2_ex_rxn (str): Reaction ID for the CO2 exchange reaction.
        Defaults to the BiGG ID 'EX_co2_e'.

    Returns:
    cue (float): CUE value
    
    """
    uptake = get_c_uptake(uptake_fluxes)
    co2_ex = get_co2_secretion(secretion_fluxes, co2_ex_rxn)
    if uptake == 0:
        cue = None
    else:
        cue = 1 - co2_ex/uptake
    
    return cue


def calculate_gge(uptake_fluxes, secretion_fluxes, co2_ex_rxn = 'EX_co2_e'):
    """Calculate the GGE by using the uptake and secretion dictionaries

    Args:
    uptake_fluxes (dict): Dictionary of carbon importing reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    secretion_fluxes (dict): Dictionary of carbon secreting reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    co2_ex_rxn (str): Reaction ID for the CO2 exchange reaction.
        Defaults to the BiGG ID 'EX_co2_e'.

    Returns:
    gge (float): GGE value
    """
    uptake = get_c_uptake(uptake_fluxes)
    co2_ex = get_co2_secretion(secretion_fluxes, co2_ex_rxn)
    org_c_ex = get_org_c_secretion(secretion_fluxes, co2_ex_rxn)
    release = co2_ex + org_c_ex
    if uptake == 0:
        gge = None
    else:
        gge = 1 - release/uptake
    
    return gge


def extract_c_fates(secretion_fluxes, uptake_fluxes = None,
                    co2_ex_rxn = 'EX_co2_e',
                    norm = False):
    """Extract the carbon atom flux going to each possible destination.
    These fluxes can the absolute value or can be normalized to the
    uptake flux.
    
    Args:
    secretion_fluxes (dict): Dictionary of carbon secreting reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    uptake_fluxes (dict): Dictionary of carbon importing reactions
        with the reaction ID and the absolute value of the carbon atom
        flux
    co2_ex_rxn (str): Reaction ID for the CO2 exchange reaction.
        Defaults to the BiGG ID 'EX_co2_e'.
    norm (bool): Do you want the flues to be normalized to the uptake
        flux? Defaults to false.

    Returns:
    
    """
    # If you want normalized fates, you must supply the uptake flux dictionary
    if norm == True and uptake_fluxes is None:
        raise ValueError('In order to calculate normalized carbon fates, ' +
                         'you must supple an uptake flux dictionary.')

    # Get the absolute values for the fate fluxes
    co2_ex = get_co2_secretion(secretion_fluxes, co2_ex_rxn)
    org_c_ex = get_org_c_secretion(secretion_fluxes, co2_ex_rxn)
    biomass_c = get_biomass_carbon()
    
    # Normalize everything to the uptake or not
    if norm == True:
        uptake = get_c_uptake(uptake_fluxes)
        if uptake == 0:
            exudation_norm = 0
            co2_ex_norm = 0
            biomass_norm = 0
        else:
            co2_ex_norm = co2_ex/uptake
            exudation_norm = org_c_ex/uptake
            biomass_norm = biomass_c/uptake
        return [co2_ex_norm, exudation_norm, biomass_norm]
    else:
        return [co2_ex, org_c_ex, biomass_c]


def extract_c_fates_from_solution(solution, c_ex_rxns, co2_ex_rxn = 'EX_co2_e', norm = True):
    # TODO: Document this function
    # Get the exchange fluxes for the current cycle
    c_ex_fluxes = {r: solution.fluxes[r] * c for r, c in c_ex_rxns.items()}
    # Use the exchange fluxes to calculate uptake, resp, and exudation
    uptake = sum([flux for rxn, flux in c_ex_fluxes.items() if flux < 0
                  and rxn != co2_ex_rxn]) # Should I count the co2_ex_rxn here?
    if c_ex_fluxes[co2_ex_rxn] < 0:
        # If the co2 flux is negative than the model is taking up CO2???
        co2_ex = 0
    else:
        co2_ex = c_ex_fluxes[co2_ex_rxn]
    exudation = abs(sum([flux for rxn, flux in c_ex_fluxes.items()
                         if flux > 0 and rxn != co2_ex_rxn]))
    # Calculate the biomass as everything that is not uptake or co2 release
    biomass = abs(uptake) - co2_ex - exudation
    # Normalize everything to the uptake or not
    if norm == True:
        if uptake == 0:
            co2_release_norm = 0
            exudation_norm = 0
            biomass_norm = 0
        else:
            co2_release_norm = co2_ex/uptake
            exudation_norm = exudation/uptake
            biomass_norm = biomass/uptake
        return [co2_release_norm, exudation_norm, biomass_norm]
    else:
        return [abs(uptake), co2_ex, exudation, biomass]
    