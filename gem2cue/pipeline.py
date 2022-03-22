import gem2cue.calculate_cue
import gem2cue.work_w_dirs
import gem2cue.plots

def pipeline(in_dir: str, out_dir: str = './results', report: bool = True, 
             boxplot: bool = True, 
             boxplot_title: str = "CUE Value for All Strains",
             boxplot_file_name: str = "CUE_boxplot",
             env_conditions_line_graphs: bool = True,
             nutrient_list = ['glc'],
             env_conditions_line_graphs_title: str = 'CUE with Variying Media and Oxygen Concentrations',
             env_conditions_line_graphs_file_name: str = 'CUE_media_O2'
            ):
    """Pipeline to do all of the analysis/make all of the figures possible
    
    Args:
        in_dir (str): Path to where all the models are saved
        out_dir (str): Path to where figures will go
        boxplot (bool): Do you want to make a boxplot of the the distribution of CUE values?
        boxplot_title (str): Title that goes on the boxplot figure
        boxplot_file_name (str): file name for the boxplot figure
        env_conditions_line_graphs (bool): Do you want to make the line graphs for the environmental conditions
        env_conditions_line_graphs_title (str): Title that goes on the line graphs figure
        env_conditions_line_graphs_file_name (str): file name for the line grpahs figure
        
    Returns:
        CUE_values (dict): Keys = model id, values = the CUE value
    """
    # List all of the files in the directory
    file_list = gem2cue.work_w_dirs.list_model_files(in_dir)

    # Read in all of the models
    model_list = gem2cue.work_w_dirs.list_cobra_models(file_list)

    # Calculate CUE for all the models
    # Make a dictionary
    CUE_values = {}
    # For each model add a key value pair of the model id and the CUE value
    for model in model_list:
        CUE_values[model.id] = gem2cue.calculate_cue.rCUE(model)

    # Get a list of just the CUE values
    data = list(CUE_values.values())

    # Start a report file
    if report:
        # Creating an HTML file
        Func = open("report.html","w")
        
        # Adding input data to the HTML file
        Func.write("<html>\n<head>\n<title> \nOutput Data in an HTML file \
                </title>\n</head> <body><h1>Welcome to <u>GeeksforGeeks</u></h1>\
                \n<h2>A <u>CS</u> Portal for Everyone</h2> \n</body></html>")
              
        # Saving the data into the HTML file
        Func.close()

    # Boxplot
    if boxplot:
        gem2cue.plots.boxplot(data, boxplot_title, out_dir, boxplot_file_name, report)

    # Environmental Line Graphs
    if env_conditions_line_graphs:
        gem2cue.plots.env_conditions_line_graphs(model_list, nutrient_list, env_conditions_line_graphs_title, out_dir, env_conditions_line_graphs_file_name, report)

    # Return values
    return(CUE_values)