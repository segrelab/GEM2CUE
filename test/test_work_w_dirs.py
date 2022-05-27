import unittest
import os

import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

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

        # Check the number of elements is as expected
        self.assertEqual(len(strain_list), 2)
        # Check the names are right
        self.assertEqual(strain_list[0].name, 'iIT341')
        self.assertEqual(strain_list[1].name, 'e_coli_core')
        # Not sure how else I can test the model
        # Check that the GC content and lengths are none
        self.assertEqual(strain_list[0].gc_content, None)
        self.assertEqual(strain_list[0].genome_length, None)
        self.assertEqual(strain_list[1].gc_content, None)
        self.assertEqual(strain_list[1].genome_length, None)

        # Add GC content and genome length
        gc_list = [0.4, 0.6]
        genome_length = [100, 200]
        strain_list = gem2cue.work_w_dirs.list_Strains(model_list,
                                                       gc_list=gc_list,
                                                       gen_length_list=genome_length)  
        self.assertEqual(len(strain_list), 2)
        self.assertEqual(strain_list[0].name, 'iIT341')
        self.assertEqual(strain_list[1].name, 'e_coli_core')
        self.assertEqual(strain_list[0].gc_content, 0.4)
        self.assertEqual(strain_list[0].genome_length, 100)
        self.assertEqual(strain_list[1].gc_content, 0.6)
        self.assertEqual(strain_list[1].genome_length, 200)

if __name__ == '__main__':
    unittest.main()