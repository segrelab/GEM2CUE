import unittest
import os
import cobra

import gem2cue.calculate_cue

class TestCalculateCUE(unittest.TestCase):
    def test_calculate_cue(self):
        """Test CUE Calculations"""
        # Read in the file to test on
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Call the function
        out_value = gem2cue.calculate_cue.rCUE(model, co2_rxn='EX_co2_e', return_sol=False)

        # Make sure that what came out is exactly what expected
        comparison_value = -9999 # FIXME: Put in a real value
        self.assertEqual(out_value, comparison_value)

if __name__ == '__main__':
    unittest.main()