from cmath import nan
import cobra
import numpy as np

def atomExchangeMetabolite(model, atom='C'):
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
    ex_atoms = {r.id: m.elements[atom] for m in model.metabolites for r in m.reactions if atom in m.elements if r.compartments == {'e'}}
    
    return ex_atoms


def rCUE(model, co2_rxn='EX_co2_e', return_sol=False):
    """ Calculate CUE, using the definition that respiration is the only waste,
    uses the formula: CUE = 1 - CO2 / Uptake C

    Args:
        model (cobra.core.Model): A model that has already been read in
        co2_rxn (string): Name of the respiration reaction in the model
        return_sol (boolean): Should the function output the FBA solution as
            well, True to return, defaults to False

    Returns:
        if return_sol = False
            outputs (int): The CUE value
        if return_sol = True
            outputs (List [int, cobra.core.Model.Solution]): The CUE and the
            last obtained solution from optimizing the model stored in a list
    """
    # Solve FBA
    sol = model.optimize()

    # Get C atoms for each exchange reaction
    ex_c_atoms = atomExchangeMetabolite(model)

    # Calculate uptake C flux
    uptake = sum([sol.get_primal_by_id(r) * ex_c_atoms[r] for r in ex_c_atoms if r != co2_rxn])
    if uptake == 0:
        cue = nan
    else:
        # Calculate CUE
        cue = 1 - abs(sol.get_primal_by_id(co2_rxn) / uptake)

    if return_sol:
        outputs = cue, sol
    else:
        outputs = cue

    return outputs


def GGE(model, return_sol=False):
    """Compute GGE using the following definition

                sum(uptake C) - sum(secretion C)
        GGE =  ----------------------------------
                        sum(uptake C)
    
    Args:
        model (cobra.core.Model): A model that has already been read in
        return_sol (boolean): Should the function output the FBA solution as
            well, True to return, defaults to False

    Returns:
        if return_sol = False
            outputs (int): The GGE value
        if return_sol = True
            outputs (List [int, cobra.core.Model.Solution]): The GGE and the
            last obtained solution from optimizing the model stored in a list
    """
    # Solve FBA
    sol = model.optimize()

    # Get C atoms for each exchange reaction
    ex_c_atoms = atomExchangeMetabolite(model)
    
    # Get C fluxes (flip signs so that uptake is positive)
    c_ex_fluxes = np.array([sol.get_primal_by_id(r) * -c for r, c in ex_c_atoms.items()])

    # Calculate GGE
    gge = c_ex_fluxes.sum() / c_ex_fluxes[c_ex_fluxes > 0].sum()

    if return_sol:
        outputs = gge, sol
    else:
        outputs = gge

    return outputs