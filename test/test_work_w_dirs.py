import unittest
import os

import gem2cue.work_w_dirs

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

class TestWorkWithDirs(unittest.TestCase):
    def test_list_model_files(self):
        """Test listing model files in a directory"""
        # Set the directory to test on
        dir_name = os.path.join(TEST_DIR, 'test_files')

        # Call the function
        file_list = gem2cue.work_w_dirs.list_model_files(dir_name)

        # Make sure that what came out is exactly what expected
        comparison_list = [TEST_DIR + '/test_files/iIT341.xml',
                           TEST_DIR + '/test_files/EC_core_flux1.xml']
        self.assertEqual(file_list, comparison_list)


    def test_list_strain_objs(self):
        """Testing listing models as Strains"""
        # Hardcode a file list to use
        file_list = [TEST_DIR + '/test_files/iIT341.xml',
                           TEST_DIR + '/test_files/EC_core_flux1.xml']

        # Read in models
        model_list = gem2cue.work_w_dirs.list_cobra_models(file_list)

        # Call the function
        strain_list = gem2cue.work_w_dirs.list_Strains(model_list)

        # Make sure that what came out is exactly what expected
        comparison_list = []
        self.assertEqual(strain_list, comparison_list)


if __name__ == '__main__':
    unittest.main()