import unittest
import os
import cobra

cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"

import gem2cue.strain

class TestStrain(unittest.TestCase):
    def test_strain_building(self):
        """Test building a Strain class object"""
        # Read in a model
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Make a Strain object
        ecoli_strain = gem2cue.strain.Strain("ecoli", model, 0.5, 100)

        # Make sure that all the Strain fields are as expected
        self.assertEqual(ecoli_strain.name, "ecoli")
        self.assertEqual(ecoli_strain.gc_content, 0.5)
        self.assertEqual(ecoli_strain.genome_length, 100)


if __name__ == '__main__':
    unittest.main()