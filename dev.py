# How to make a boxplot with two boxes on it

import matplotlib.pyplot as plt

my_dict = {'ABC': [34.54, 34.345, 34.761], 'DEF': [34.541, 34.748, 34.482]}

fig, ax = plt.subplots()
ax.boxplot(my_dict.values())
ax.set_xticklabels(my_dict.keys())

plt.savefig("test.png")