import unittest
import os
import shutil
import cobra

import gem2cue.plots
import gem2cue.work_w_dirs

import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(TEST_DIR, 'test_results')

class TestPlots(unittest.TestCase):
    def test_sankey_plot(self):
        "Test Sankey Plot function"
        # Make sure there is not an output direcotory with that name already
        if os.path.exists(OUT_DIR):
            raise ValueError(f'There is already a directory named {OUT_DIR}')
        # Then make it
        os.makedirs(OUT_DIR)

        # Make plot for ecoli
        model_file = os.path.join(TEST_DIR, 'test_files', 'EC_core_flux1.xml')
        model = cobra.io.read_sbml_model(model_file)
        gem2cue.plots.sankey_plot(model, out_dir = OUT_DIR)

        # Make plot for other test model
        model_file = os.path.join(TEST_DIR, 'test_files', 'iIT341.xml')
        model = cobra.io.read_sbml_model(model_file)
        gem2cue.plots.sankey_plot(model, biomass_rxn = 'BIOMASS_HP_published',
                                  out_dir = OUT_DIR)

        # Look at the number of figures check the number is as expected
        n_figs = len([name for name in os.listdir(OUT_DIR)])
        self.assertEqual(n_figs, 2)


    def test_gc_len_plot(self):
        "Test plots of CUE vs GC content and genome length"
        # Make sure there is not an output direcotory with that name already
        if os.path.exists(OUT_DIR):
            raise ValueError(f'There is already a directory named {OUT_DIR}')
        # Then make it
        os.makedirs(OUT_DIR)

        # Read in the models and make a list of strains
        file_1 = os.path.join(TEST_DIR, 'test_files', 'EC_core_flux1.xml')
        model_1 = cobra.io.read_sbml_model(file_1)

        file_2 = os.path.join(TEST_DIR, 'test_files', 'iIT341.xml')
        model_2 = cobra.io.read_sbml_model(file_2)
        
        strain_list = gem2cue.work_w_dirs.list_Strains([model_1, model_2])

        # Make GC plot
        gem2cue.plots.gc_content(strain_list,
                                 definition = gem2cue.calculate_cue.rCUE,
                                 title = 'CUE vs GC content', 
                                 out_dir = OUT_DIR,
                                 file_name = 'gc_plot',
                                 report = False)

        # Make genome length plot
        gem2cue.plots.genome_length(strain_list,
                                    definition = gem2cue.calculate_cue.rCUE,
                                    title = 'CUE vs GC content', 
                                    out_dir = OUT_DIR,
                                    file_name = 'gen_len_plot',
                                    report = False)

        # Look at the number of figures check the number is as expected
        n_figs = len([name for name in os.listdir(OUT_DIR)])
        self.assertEqual(n_figs, 2)

        
    def tearDown(self):
        # Delete the output folder
        shutil.rmtree(OUT_DIR)


if __name__ == '__main__':
    unittest.main()