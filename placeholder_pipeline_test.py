# I am not sure why VSCode won't show me the pipeline test
# Recapitualted the code here for easy running

import gem2cue.pipeline

# nutrients = ['glc', 'accoa_c', 'fum_c']
gem2cue.pipeline.pipeline('./test/test_files', './test/test_results', boxplot=False)