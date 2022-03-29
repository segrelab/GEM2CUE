# Basically re-do the pipeline with strain objects
import os
import cobra
import matplotlib.pyplot as plt

import gem2cue.strain
import gem2cue.calculate_cue

# Read in the models
test_dir = '/home/helen/Documents/PhD/Segre Lab/GEM2CUE/test'
ecoli_model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))
ilT341_model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'iIT341.xml'))

# Make strain objects for my two test files
ecoli_strain = gem2cue.strain.Strain("ecoli", ecoli_model, 0.5, 4641652)
ilT341_strain = gem2cue.strain.Strain("ilT341", ilT341_model, 0.3, 1667867)
strain_list = [ecoli_strain, ilT341_strain]

# Create lists of the CUE and the metadata
gc_list = []
length_list = []
CUE_list = []
for strain in strain_list:
    gc_list.append(strain.gc_content)
    length_list.append(strain.genome_length)
    CUE_list.append(gem2cue.calculate_cue.rCUE(strain.model))

# Plot against GC content
plt.plot(gc_list, CUE_list)
plt.xlabel('GC Content (%)')
plt.ylabel('CUE')
plt.savefig("test_CUE_vs_GC.png")
plt.close()

# Plot against genome length
plt.plot(length_list, CUE_list)
plt.xlabel('Genome Length (bp)')
plt.ylabel('CUE')
plt.savefig("test_CUE_vs_length.png")
plt.close()