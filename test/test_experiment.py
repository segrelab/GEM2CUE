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

        # Make a Media object
        medium = gem2cue.utils.Media({'EX_co2_e': 1000.0,
            'EX_glc__D_e': 10.0,
            'EX_h_e': 1000.0,
            'EX_h2o_e': 1000.0,
            'EX_nh4_e': 1000.0,
            'EX_pi_e': 1000.0})

        # Make an Experiment object
        ecoli_exp = gem2cue.utils.Experiment(ecoli_strain, medium)

        # Make sure that all the Experiment attributes are as expected
        self.assertEqual(ecoli_exp.strain.name, "ecoli")
        self.assertEqual(ecoli_exp.strain.metadata['gc_content'], 0.5)
        self.assertEqual(ecoli_exp.strain.metadata['genome_length'], 100)

        self.assertEqual(ecoli_exp.media.media['EX_co2_e'], 1000.0)
        self.assertEqual(ecoli_exp.media.media['EX_glc__D_e'], 10.0)
        self.assertEqual(ecoli_exp.media.media['EX_h_e'], 1000.0)

        # Run the "run" method
        ecoli_exp.run()

        # Check that the solution is as expected
        self.assertEqual(ecoli_exp.solution.objective_value, 0.8739215069684279)

        # Run the CUE method
        ecoli_exp.CUE()

        # Check that the result is as expected
        self.assertEqual(ecoli_exp.cue, 0.6198361114965819)


    def test_media_change(self):
        """Test that supplying a media to an experiment will overwrite the 
        COBRApy medium"""
        # Read in a model
        test_dir = os.path.dirname(os.path.realpath(__file__))
        model = cobra.io.read_sbml_model(os.path.join(test_dir, 'test_files', 'EC_core_flux1.xml'))

        # Make a Strain object
        ecoli_strain = gem2cue.utils.Strain("ecoli", model, {'gc_content' :0.5, 'genome_length': 100})

        # Define a media object
        medium = gem2cue.utils.Media({'EX_co2_e': 1000.0,
            'EX_glc__D_e': 10.0,
            'EX_h_e': 1000.0,
            'EX_h2o_e': 1000.0,
            'EX_nh4_e': 1000.0,
            'EX_pi_e': 1000.0})

        # Define an experiment
        ecoli_exp = gem2cue.utils.Experiment(ecoli_strain, medium)

        # Check that the COBRApy medium still hasn't changed
        self.assertEqual(ecoli_exp.strain.model.medium, 
                        {'EX_co2_e': 1000.0,
                        'EX_glc__D_e': 10.0,
                        'EX_h_e': 1000.0,
                        'EX_h2o_e': 1000.0,
                        'EX_nh4_e': 1000.0,
                        'EX_o2_e': 1000.0,
                        'EX_pi_e': 1000.0})

        # Run FBA
        ecoli_exp.run()

        # Check that the medium changed to what I specified
        self.assertEqual(ecoli_exp.strain.model.medium, 
                        {'EX_co2_e': 1000.0,
                        'EX_glc__D_e': 10.0,
                        'EX_h_e': 1000.0,
                        'EX_h2o_e': 1000.0,
                        'EX_nh4_e': 1000.0,
                        'EX_pi_e': 1000.0})



if __name__ == '__main__':
    unittest.main()