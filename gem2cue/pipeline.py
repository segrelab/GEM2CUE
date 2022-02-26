import gem2cue.calculate_cue
import gem2cue.work_w_dirs
import gem2cue.plots

def pipeline(in_dir: str, out_dir:str = './results', boxplot=True, 
             boxplot_title = "CUE Value for All Strains",
             boxplot_file_name = "CUE_boxplot"):
    # List all of the files in the directory
    file_list = gem2cue.work_w_dirs.list_model_files(in_dir)

    # Read in all of the models
    model_list = gem2cue.work_w_dirs.list_cobra_models(file_list)

    # Calculate CUE for all the models
    # Make a dictionary
    CUE_values = {}
    # For each model add a key value pair of the model name and the CUE value
    for model in model_list:
        CUE_values[model.name] = gem2cue.calculate_cue.rCUE(model)

    # Get a list of just the CUE values
    data = list(CUE_values.values)

    # Boxplot
    gem2cue.plots.boxplot(data, boxplot_title, out_dir, boxplot_file_name)