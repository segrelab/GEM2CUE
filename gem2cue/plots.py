import matplotlib as plt

def boxplot(data, title, out_dir, file_name):
    """Take a list of CUE values, plot them as a boxplot, and save to the 
    specified output directory

    Args:
        data ([float]): List of CUE values to plot
        title (str): Title of the boxplot
        out_dir (str): Location to save the boxplot
        file_name (str): Name of the file to save

    Returns:
        Nothing, saves figure file to the output directory
    """
    fig = plt.figure(figsize =(10, 7))
 
    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])

    # Adding title
    ax.set_title(title)
    
    # Creating plot
    ax.boxplot(data)

    # Axis Labels
    ax.set_ylabel('CUE')

    # If the output directory does not exsit, make it
    import os
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save
    plt.savefig(out_dir + "/" + file_name + ".png")
    plt.close()