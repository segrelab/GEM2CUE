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
        comparison_value = 0.6198361114965837
        self.assertEqual(out_value, comparison_value)


    def test_atomExchangeMetabolite(self):
        """Test the directory function"""
        # Read in the file to test on
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Call the function
        out_value = gem2cue.calculate_cue.atomExchangeMetabolite(model)

        # Make sure that what came out is exactly what expected
        comparison_value = {'EX_ac_e': 2,
                            'EX_acald_e': 2,
                            'EX_akg_e': 5,
                            'EX_co2_e': 1,
                            'EX_etoh_e': 2,
                            'EX_for_e': 1,
                            'EX_fru_e': 6,
                            'EX_fum_e': 4,
                            'EX_glc__D_e': 6,
                            'EX_gln__L_e': 5,
                            'EX_glu__L_e': 5,
                            'EX_lac__D_e': 3,
                            'EX_mal__L_e': 4,
                            'EX_pyr_e': 3,
                            'EX_succ_e': 4}
        self.assertEqual(out_value, comparison_value)

    
    def test_calculate_gge(self):
        """Test GGE Calculations"""
        # Read in the file to test on
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Call the function
        out_value = gem2cue.calculate_cue.GGE(model, co2_rxn='EX_co2_e', return_sol=False)

        # Make sure that what came out is exactly what expected
        comparison_value = -9999
        self.assertEqual(out_value, comparison_value)


if __name__ == '__main__':
    unittest.main()