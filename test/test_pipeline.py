import unittest
import os
import shutil

import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

import gem2cue.pipeline

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
OUT_DIR = os.path.join(TEST_DIR, 'test_results')

class TestPipeline(unittest.TestCase):
    def test_pipeline(self):
        # Set the directory to pass
        # Set the directory to test on
        dir_name = os.path.join(TEST_DIR, 'test_files')

        # Make sure there is not an output direcotory with that name already
        if os.path.exists(OUT_DIR):
            raise ValueError(f'There is already a directory named {OUT_DIR}')

        # Run the pipeline
        gem2cue.pipeline.pipeline(dir_name, OUT_DIR)

        # Look at the number of figures (and report), check the number is as expected
        n_figs = len([name for name in os.listdir(OUT_DIR)])
        self.assertEqual(n_figs, 5)

        
    def tearDown(self):
        # Delete the output folder
        shutil.rmtree(OUT_DIR)


if __name__ == '__main__':
    unittest.main()