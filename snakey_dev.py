import matplotlib.pylab as plt
from matplotlib.sankey import Sankey
import cobra

import gem2cue.calculate_cue

################################################################################
# Example 1 -- Mostly defaults
#
# This demonstrates how to create a simple diagram by implicitly calling the
# Sankey.add() method and by appending finish() to the call to the class.

# Sankey(flows=[0.25, 0.15, 0.60, -0.20, -0.15, -0.05, -0.50, -0.10],
#        labels=['', '', '', 'First', 'Second', 'Third', 'Fourth', 'Fifth'],
#        orientations=[-1, 1, 0, 1, 1, 1, 0, -1]).finish()
# plt.title("The default settings produce a diagram like this.")

# plt.savefig("sankey_ex_01.png")
################################################################################
# Made up FBA results
# Sankey(flows = [1, -0.33, -0.33, -0.33],
#        labels = ['U = A', 'R', 'EX', 'G'],
#        orientations = [0, 1, 1, 0],
#        pathlengths = [0.6, 0.25, 0.25, 0.4]).finish()
# plt.title("CUE definition")
# plt.savefig("sankey_ex_02.png")
################################################################################
# Real FBA results
# Read in the model and solve it
model_file = 'test/test_files/EC_core_flux1.xml'
model = cobra.io.read_sbml_model(model_file)
sol = model.optimize()

# Reactions from the model you need to know the name of
co2_rxn = 'EX_co2_e'
biomass_rxn = 'BIOMASS_Ecoli_core_w_GAM'

# Look up the exchange reactions
ex_c_atoms = gem2cue.calculate_cue.atomExchangeMetabolite(model)

# Find the uptake flux
uptake_rxns = sol.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x<0]).index # Negative flux for uptake
uptake = abs(sum([sol.get_primal_by_id(r) * ex_c_atoms[r] for r in uptake_rxns if r != co2_rxn])) # Doesn't really matter if we say no CO2 rxn or not, it shouldn't be a negative flux

# Find the respiration flux
resp = -1 * sol.get_primal_by_id(co2_rxn) # Don't multiply because we know that the number of c atoms is 1 for CO2

# Find the exudation flux
exudation_rxns = sol.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x>0]).index # Postive flux for exudation
ex = abs(sum([sol.get_primal_by_id(r) * ex_c_atoms[r] for r in exudation_rxns if r != co2_rxn])) # Do really need to say no CO2 reaction

# Find the biomass flux
# Make a dictionary of the fluxes for each metabolite in the biomass reaction
biomass_dict = model.reactions.BIOMASS_Ecoli_core_w_GAM.metabolites
# Find the number of C atoms that the biomass reaction takes in
biom_in = abs(sum([biomass_dict[m] * sol.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] < 0 if 'C' in m.elements]))
# Find the number of C atoms that the biomass reaction puts out
biom_out = abs(sum([biomass_dict[m] * sol.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] > 0 if 'C' in m.elements]))
# Find the difference
# Is this what actually goes to biomass?
biom = -1 * (biom_in - biom_out)

# Make the ex a very small negative number, just to try to fix the graph shape
# ex = -0.01

# Make the plot
Sankey(flows = [uptake, resp, ex, biom],
       labels = ['U = A', 'R', 'EX', 'G'],
       orientations = [0, 1, 1, 0],
       pathlengths = [50, 10, 10, 20]).finish()
plt.title("FBA results for E. coli core model")
plt.savefig("sankey_ex_03.png")
