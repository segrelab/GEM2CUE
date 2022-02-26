import unittest
import os

import gem2cue.work_w_dirs

class TestWorkWithDirs(unittest.TestCase):
    def test_list_model_files(self):
        """Test listing model files in a directory"""
        # Set the directory to test on
        test_dir = os.path.dirname(os.path.realpath(__file__))
        dir_name = os.path.join(test_dir, 'test_files')

        # Call the function
        file_list = gem2cue.work_w_dirs.list_model_files(dir_name)

        # Make sure that what came out is exactly what expected
        comparison_list = [test_dir + '/test_files/iLJ478.xml',
                           test_dir + '/test_files/iIT341.xml',
                           test_dir + '/test_files/EC_core_flux1.xml']
        self.assertEqual(file_list, comparison_list)


if __name__ == '__main__':
    unittest.main()