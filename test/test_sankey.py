import unittest
import os
import shutil
import cobra

import gem2cue.plots

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(TEST_DIR, 'test_results')

class TestPipeline(unittest.TestCase):
    def test_sankey_plot(self):
        # Read in a model
        model_file = os.path.join(TEST_DIR, 'test_files', 'EC_core_flux1.xml')
        model = cobra.io.read_sbml_model(model_file)

        # Make sure there is not an output direcotory with that name already
        if os.path.exists(OUT_DIR):
            raise ValueError(f'There is already a directory named {OUT_DIR}')

        # Run the pipeline
        gem2cue.plots.sankey_plot(model, out_dir = OUT_DIR)

        # Look at the number of figures check the number is as expected
        n_figs = len([name for name in os.listdir(OUT_DIR)])
        self.assertEqual(n_figs, 1)

        
    def tearDown(self):
        # Delete the output folder
        shutil.rmtree(OUT_DIR)


if __name__ == '__main__':
    unittest.main()