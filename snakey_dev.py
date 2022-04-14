import matplotlib.pylab as plt
from matplotlib.sankey import Sankey

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
Sankey(flows = [1, -0.33, -0.33, -0.33],
       labels = ['U = A', 'R', 'EX', 'G'],
       orientations = [0, 1, 1, 0],
       pathlengths = [0.6, 0.25, 0.25, 0.4]).finish()
plt.title("CUE definition")
plt.savefig("sankey_ex_02.png")


