import unittest
import os
import shutil

import gem2cue.pipeline

class TestPipeline(unittest.TestCase):
    def test_pipeline(self):
        # Set the directory to pass
        # Set the directory to test on
        test_dir = os.path.dirname(os.path.realpath(__file__))
        dir_name = os.path.join(test_dir, 'test_files')

        # Make sure there is not an output direcotory with that name already
        out_dir = os.path.join(test_dir, 'test_results')
        if os.path.exists(out_dir):
            raise ValueError(f'There is already a directory named {out_dir}')

        # Run the pipeline
        gem2cue.pipeline.pipeline(dir_name, out_dir)

        # Look at the number of figures, check the number is as expected
        n_figs = len([name for name in os.listdir(out_dir)])
        self.assertEqual(n_figs, 1)

        # Delete the output folder
        shutil.rmtree(out_dir)


if __name__ == '__main__':
    unittest.main()