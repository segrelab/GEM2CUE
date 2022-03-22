import glob
import cobra

def list_model_files(dir: str, type: str = 'xml'):
    """Take a directory and list all of the model files in it
    
    Args:
        dir (str): name of directory to search
        type (str): file type (extension) to filter for
        
    Returns:
        file_list ([str]): List of all the model names
    """
    
    file_list = [f for f in glob.glob(dir + "/*." + type)]

    return(file_list)

def list_cobra_models(file_list):
    """Take a list of files and read them in as cobra models
    
    Args:
        file_list ([str]): List of all the model names
        
    Returns:
        model_list ([cobra.Model]): List of all the model objects

    Assumptions/Limitations:
        Currently only works for SMBL files
    """
    # TODO:
    # Check what the file endings are
    # Pick the right function to use
    
    model_list = [cobra.io.read_sbml_model(f) for f in file_list]

    return(model_list)