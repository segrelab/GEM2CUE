import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

import gem2cue.calculate_cue

def boxplot(data, title, out_dir, file_name, report):
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
    # Make the figure, create axes instance
    fig, ax = plt.subplots()

    # Adding title
    ax.set_title(title)
    
    # Creating plot
    ax.boxplot(data)

    # Axis Labels
    ax.set_ylabel('CUE')
    ax.axes.xaxis.set_visible(False)

    # If the output directory does not exsit, make it
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save
    plt.savefig(out_dir + "/" + file_name + ".png")
    plt.close()
    
    # Add to the report
    if report:
        pass


def env_conditions_line_graphs(model_list, nutrient_list, title, out_dir, file_name, report):
    """
    
    Args:

    Returns:

    """
    # Make a header for the results in the report
    if report:
        pass

    # Make a set of plots for each individiual model
    for model in model_list:
        # Get the name of the model, if there is no name- call it no name
        # TODO: Throw a warning because if you have multiple models with no name you will overwrite results
        if model.name == '':
            name = 'Unknown Model'
        else:
            name = model.name

        # Save CUE for each nutrient and o2 concentration pair in a data frame
        data = []
        for nutrient in nutrient_list:
            for nutrient_conc in range(10, 21):
                # Update the glucose in the medium
                medium = model.medium
                # Need to look up the reaction name to do nutrients other than glucose
                medium['EX_glc__D_e'] = nutrient_conc
                for o2 in np.linspace(0, 25, 5):
                    # Update the oxygen in the medium
                    medium['EX_o2_e'] = o2
                    # Set the model to use the updated medium
                    model.medium = medium
                    # Calculate CUE
                    cue, sol = gem2cue.calculate_cue.rCUE(model, return_sol = True)

                    # Save to the data frame
                    row = {'nutrient': nutrient, 'nutrient_conc': nutrient_conc,
                           'o2': o2, 'biomass': sol.objective_value,
                           'co2': sol.fluxes.EX_co2_e, 'CUE': cue}
                    data.append(row)

        # Convert results to a data fram
        df = pd.DataFrame(data)

        # Plot
        g = sns.relplot(x='nutrient_conc', y='CUE', hue='nutrient', col='o2', data=df, kind='line', marker='o', height=3)

        # Add overall title to replot
        g.fig.suptitle(title + " for " + name)

        # Adjust subplots so that titles don't overlap
        g.fig.subplots_adjust(top = 0.85)

        # Set the x axis label
        g.set_xlabels("Nutrient Concentration\n (mol / [1 gDW * 24 h])")
        
        # If the output directory does not exsit, make it
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Save
        plt.savefig(out_dir + "/" + file_name + "_" + name + ".png")
        plt.close()

        # Add to the report
        if report:
            pass