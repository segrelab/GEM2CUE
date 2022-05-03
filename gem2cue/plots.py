import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import numpy as np
import pandas as pd
import seaborn as sns
import os
import cobra

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
        Func = open(os.path.join(out_dir, "report.html"), "a")
        Func.write('\n<img src="' + file_name + '.png">')
        Func.close()


def env_conditions_line_graphs(model_list, nutrient_list, definition, title, out_dir, file_name, report):
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
                    cue, sol = definition(model, return_sol = True)

                    # Save to the data frame
                    row = {'nutrient': nutrient, 'nutrient_conc': nutrient_conc,
                           'o2': o2, 'biomass': sol.objective_value,
                           'co2': sol.fluxes.EX_co2_e, 'CUE': cue}
                    data.append(row)

        # Convert results to a data fram
        df = pd.DataFrame(data)

        # Plot
        # FIXME: Make the y axis name based on the definition you are using
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
            Func = open(os.path.join(out_dir, "report.html"), "a")
            Func.write('\n<img src="' + file_name + '_' + name + '.png">')
            Func.close()


def sankey_plot(model: cobra.core.Model,
           co2_rxn: str = 'EX_co2_e',
           biomass_rxn: str = 'BIOMASS_Ecoli_core_w_GAM',
           title_stem: str = 'Carbon Flow for ',
           out_dir: str = '.',
           file_stem: str = 'sankey_'):
    """ Create and save a sankey diagram for the CUE of a specific model
    
    Args:
    model (cobra.core.Model): Model to be solved
    co2_rxn (str): Name of the respiration reaction in the model (defaults to 
        'EX_co2_e')
    biomass_rxn (str): Name of the biomass reaction in the model (defaults to
        'BIOMASS_Ecoli_core_w_GAM')
    title_stem (str): First part of the title on the plot, will be catted with
        the model name (defaults to 'Carbon Flow for ')
    out_dir (str): Output path (defaults to '.')
    file_stem (str): First part of the file name, will be catted with the model
        name (defaults to 'sankey_')
    
    Returns:
    Nothing, figure is automatically saved
    """
    # Solve the model
    sol = model.optimize()

    # Look up the exchange reactions
    ex_c_atoms = gem2cue.calculate_cue.atomExchangeMetabolite(model)

    # Find the uptake flux
    uptake_rxns = sol.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x<0]).index # Negative flux for uptake
    uptake = abs(sum([sol.get_primal_by_id(r) * ex_c_atoms[r] for r in uptake_rxns if r != co2_rxn])) # Doesn't really matter if we say no CO2 rxn or not, it shouldn't be a negative flux

    # Find the respiration flux
    resp = -1 * sol.get_primal_by_id(co2_rxn) # Don't multiply because we know that the number of c atoms is 1 for CO2

    # Find the exudation flux
    exudation_rxns = sol.fluxes[ex_c_atoms.keys()].pipe(lambda x: x[x>0]).index # Postive flux for exudation
    ex = -1 * abs(sum([sol.get_primal_by_id(r) * ex_c_atoms[r] for r in exudation_rxns if r != co2_rxn])) # Do really need to say no CO2 reaction

    # Find the biomass flux
    # Make a dictionary of the fluxes for each metabolite in the biomass reaction
    biomass_dict = getattr(model.reactions, biomass_rxn).metabolites
    # Find the number of C atoms that the biomass reaction takes in
    biom_in = abs(sum([biomass_dict[m] * sol.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] < 0 if 'C' in m.elements]))
    # Find the number of C atoms that the biomass reaction puts out
    biom_out = abs(sum([biomass_dict[m] * sol.get_primal_by_id(biomass_rxn) * m.elements['C'] for m in biomass_dict if biomass_dict[m] > 0 if 'C' in m.elements]))
    # Find the difference
    # Is this what actually goes to biomass?
    biom = -1 * (biom_in - biom_out)


    # Make the plot
    Sankey(flows = [uptake, resp, ex, biom],
        labels = ['U = A', 'R', 'EX', 'G'],
        orientations = [0, 1, 1, 0],
        pathlengths = [50, 10, 10, 20]).finish()
    # Make sure the title will have a model name in it
    if model.name == '':
        print_name = model.id
    else:
        print_name = model.name
    plt.title(title_stem + print_name)
    plt.savefig(os.path.join(out_dir, file_stem + model.id + ".png"))