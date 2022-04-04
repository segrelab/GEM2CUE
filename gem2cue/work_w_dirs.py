import glob
import cobra

import gem2cue.strain

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


def list_Strains(model_list, name_list = None, gc_list = None, gen_length_list = None):
    """Take a list of models and make them Strain objects, if any metadata is
    provided, include that in the Strain
    
    Args:
        model_list ([cobra.Model]): List of all the model objects
        gc_list ([int]): GC content for the strains
        gen_length_list ([int]): Genome lengths for the strains
        
    Returns:
        strain_list ([Strain]): List of all the strain objects
    """
    # TODO: Ensure that all the input lists are the same length
    strain_list = []
    for idx in range(len(model_list)):
        model = model_list[idx]
        # Collect the metadata
        if name_list is None:
            name = model.id
        if gc_list is None:
            gc = None
        else:
            gc = gc_list[idx]
        if gen_length_list is None:
            length = None
        else:
            length = gen_length_list[idx]

        # Make the strain object
        strain_obj = gem2cue.strain.Strain(name,
                                           model,
                                           gc,
                                           length)
        strain_list.append(strain_obj)
    
    return(strain_list)