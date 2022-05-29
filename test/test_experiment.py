from operator import mod
import unittest
import os
import cobra

import gem2cue.utils

class TestExperiment(unittest.TestCase):
    def test_experiment(self):
        """Test building an Experiment class object and running its methods"""
        # Read in a model
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Make a Strain object
        ecoli_strain = gem2cue.utils.Strain("ecoli", model, {'gc_content' :0.5, 'genome_length': 100})

        # Make an Experiment object
        ecoli_exp = gem2cue.utils.Experiment(ecoli_strain)

        # Make sure that all the Experiment attributes are as expected
        self.assertEqual(ecoli_exp.strain.name, "ecoli")
        self.assertEqual(ecoli_exp.strain.metadata['gc_content'], 0.5)
        self.assertEqual(ecoli_exp.strain.metadata['genome_length'], 100)

        # Run the "run" method
        ecoli_exp.run()

        # Check that the solution is as expected
        self.assertEqual(ecoli_exp.solution.objective_value, 0.8739215069684279)

        # Run the CUE method
        ecoli_exp.CUE()

        # Check that the result is as expected
        self.assertEqual(ecoli_exp.cue, 0.6198361114965819)

if __name__ == '__main__':
    unittest.main()