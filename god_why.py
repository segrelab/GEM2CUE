import os
import cobra

import gem2cue.utils

# Make a strain
test_strain = gem2cue.utils.Strain('test',
	cobra.io.read_sbml_model(os.path.join('IAMM_models', 'faa_from_ncbi', 'GCF_000011965.2.xml')),
	None,
	None)

# Make an experiment
glc_10_exp = gem2cue.utils.Experiment(test_strain)

# Edit the medium to have lower glucose
medium = glc_10_exp.strain.model.medium
medium['EX_glc__D_e'] = 10
glc_10_exp.strain.model.medium = medium

# Run
glc_10_exp.CUE(ex_nomenclature = {'C_e'})
print(glc_10_exp.cue)

# Make an experiment
no_glc_exp = gem2cue.utils.Experiment(test_strain)

# Edit the medium to have lower glucose
medium = no_glc_exp.strain.model.medium
medium['EX_glc__D_e'] = 0
no_glc_exp.strain.model.medium = medium

# Run
no_glc_exp.CUE(ex_nomenclature = {'C_e'})
print(no_glc_exp.cue)