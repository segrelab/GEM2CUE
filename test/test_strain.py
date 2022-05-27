import unittest
import os
import cobra

import gem2cue.utils

class TestStrain(unittest.TestCase):
    def test_strain(self):
        "Test building a Strain class object and changing the media"
        # Read in a model
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Make a Strain object
        ecoli_strain = gem2cue.utils.Strain("ecoli", model, {'gc_content' :0.5, 'genome_length': 100})

        # Make sure that all the Strain fields are as expected
        self.assertEqual(ecoli_strain.name, "ecoli")
        self.assertEqual(ecoli_strain.metadata['gc_content'], 0.5)
        self.assertEqual(ecoli_strain.metadata['genome_length'], 100)

        # Make sure that the original medium is the unconstrained
        self.assertEqual(len(ecoli_strain.model.medium), 7)
        self.assertEqual(list(ecoli_strain.model.medium.values()),
            list(ecoli_strain.model.medium.values()))

        # Make a new medium
        new_med = gem2cue.utils.Media({'EX_co2_e': 1000.0, 'EX_ac_e': 10.0, 'EX_h_e': 1000.0, 'EX_h2o_e': 1000.0, 'EX_nh4_e': 1000.0, 'EX_pi_e': 1000.0})

        # Change the medium for the strain
        ecoli_strain.update_medium(new_med)

        # Make sure that the new medium is as expected
        self.assertEqual(ecoli_strain.model.medium,
            {'EX_ac_e': 10.0, 'EX_co2_e': 1000.0, 'EX_h_e': 1000.0, 'EX_h2o_e': 1000.0, 'EX_nh4_e': 1000.0, 'EX_pi_e': 1000.0})


if __name__ == '__main__':
    unittest.main()