#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot altitude profiles for a variety of run outputs into PDF files
"""

__author__ = "Daniel Lo, Roger Yelle"

#%% ----------------------------- Packages Import -----------------------------

import sys
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfWriter, PdfReader
import saguaro_read as sgr


#%% ---------------------- Class and Method Definitions -----------------------

# quantities for each plot
class plot_variables:
    def __init__(self, variables, xlabel):
        self.variables = variables
        self.xlabel = xlabel
    def append_variable(self, x):
        self.variables.append(x)


# returns a string of the reaction
def reaction_title(r1, r2, p1, p2, p3, p4):
    title = r1 + ' + ' + r2 + ' → ' + p1
    if p2 != '': title = title + ' + ' + p2
    if p3 != '': title = title + ' + ' + p3
    if p4 != '': title = title + ' + ' + p4
    return title


# get and sort reactions by column rates
def get_sorted_reactions(run_name):
    # read reaction rates
    crates = sgr.rates('../runs/'+run_name+'/output/chemrates.out')
    prates = sgr.rates('../runs/'+run_name+'/output/photorates.out')
    erates = sgr.rates('../runs/'+run_name+'/output/elerates.out')
    alt = crates['Altitude']
    n_alt = len(alt)
    
    reactions = pd.concat(
        [crates['Reaction'], prates['Reaction'], erates['Reaction']],
        ignore_index=True)
    n_rct = len(reactions)
    
    rates = np.concatenate((crates['Rates'], prates['Rates'], erates['Rates']),
                           axis=1)
    
    # create table of column integrated rates for filtering and sorting
    R = 3390 # Mars radius
    rz = R + alt
    rz2 = np.power(rz/R, 2)
    colrates = np.zeros(n_rct)
    for i_r in range(0, n_rct):
       sm = 0.
       for i_z in range(1, n_alt):
           sm = sm + ((rz[i_z] - rz[i_z-1]) 
                      * (rates[i_z, i_r] * rz2[i_z] 
                         + rates[i_z-1, i_r] * rz2[i_z-1]) / 2)
       colrates[i_r] = sm
    reactions['Column Rate'] = colrates
    reactions_grouped = reactions.groupby(['Reactant1', 'Reactant2', 
                                           'Product1', 'Product2',
                                           'Product3', 'Product4'],
                                          sort=False)
    sort_order = reactions_grouped['Column Rate'].sum().sort_values(
        ascending=False)
    return alt, rates, reactions_grouped, sort_order


# Make 4x4 altitude-profile plots and save into multi-page PDF.
# plot_sets: sets of quantities to be included in each plot
# data: Pandas DataFrame with altitude profile values for each heading
def plot_profiles(save_location, alt, plot_sets, data, plot_title=''):
    with PdfPages(save_location) as pdf:
        cmap = plt.cm.Set1
        i_plots = 0
        n_plots = len(plot_sets)
        n_pages = int(np.ceil(n_plots/4))
        for i_page in range(n_pages):
            fig = plt.figure(figsize=[8,8], layout='compressed')
            fig.suptitle(plot_title, fontsize=14)
            for i_subplot in range(1,5):
                variables = plot_sets[i_plots].variables
                ax = fig.add_subplot(2, 2, i_subplot)
                for i in range(len(variables)):
                    color = cmap(i / (cmap.N - 1))
                    ax.semilogx(data[variables[i]], alt, 
                                color = color, label = variables[i])
                    # ax.text(0.98, 0.98-i*0.03, 
                    #         variables[i],
                    #         horizontalalignment = 'right', 
                    #         verticalalignment = 'top',
                    #         transform = ax.transAxes,
                    #         color = color, fontsize = '6')

                ax.legend(fontsize=6, frameon=False, 
                          labelcolor='linecolor', 
                          handlelength=0, handletextpad=0)
                xmax = log10_ceil(np.max(data[variables]))*10
                xmin = max([
                    log10_floor(np.min(data[variables])),
                    log10_floor(np.min(np.max(data[variables], axis=0))) / 1.E8
                    ])
                ymin = 20 * np.floor(min(alt) / 20) 
                ymax = 50 * np.ceil(max(alt) / 50)
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_xlabel(plot_sets[i_plots].xlabel)
                if (i_subplot == 1) or (i_subplot == 3):
                    ax.set_ylabel('Altitude (km)')
                    
                i_plots += 1  
                if i_plots == n_plots:
                    break
            pdf.savefig()
            plt.close()
            if i_plots == n_plots:
                break
        print('Plots saved to ' + save_location + '.')


# ceiling in log10-space
def log10_ceil(x):
    return 10 ** np.ceil(np.log10(x))

# floor in log10-space
def log10_floor(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return 10 ** np.floor(np.log10(x))

# split positive and negative values into upward and downward series
def split_up_down(x):
    up = x.copy()
    up[up < 0] = np.nan
    down = -x.copy()
    down[down < 0] = np.nan
    return up, down


#%% -------------------------------- Run Code ---------------------------------

# Get user input on data and plots
try: 
    p = int(input('Select figure to plot:\n'
          +'    [1] Densities\n'
          +'    [2] Mixing ratios\n'
          +'    [3] Single species\n'
          +'    [4] All reactions rates\n'
          +'>>> '))
    if (p<1) or (p>4):
        raise ValueError
except: 
    print('Invalid input! Exiting...')
    sys.exit()

run_name = input('Enter run name > ').strip()
if not os.path.isdir('../runs/' + run_name):
    print(f"Error: Run '{run_name}' does not exist.")
    sys.exit()

plot_path = '../runs/' + run_name + '/plots/'


if p == 1: # plot densities
    # Read in atm1D.out
    atm = sgr.atm('../runs/' + run_name + '/output/atm1D.out')
    alt = 1.E-5 * atm['Profiles']['Altitude'] # convert cm to km
    
    # Plot neutrals and save to temporary file
    plot_sets = []
    t_plot_set = plot_variables(variables=[], xlabel=r'Density (cm$^{-3}$)')
    plot_set_index = 0
    for i in range(len(atm['Species'])):
        if 'P' in atm['Species'][i]: # stop when encounter ion
            if t_plot_set.variables: plot_sets.append(t_plot_set)
            ion_index = i
            break
        t_plot_set.append_variable(atm['Species'][i])
        plot_set_index += 1
        
        # reach 5 variables for plot or run out of variables
        if (plot_set_index == 5) or (i == len(atm['Species'])-1):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel=r'Density (cm$^{-3}$)')
            plot_set_index = 0
    
    pdf_neutral_name = 'density_temp_n.pdf'
    plot_profiles(plot_path + pdf_neutral_name, 
                  alt, plot_sets, atm['Profiles'], 
                  plot_title = 'Number Density (Neutrals)')

    # Plot electrons and ions and save to temporary file
    # move electrons to front
    plot_sets = []
    t_plot_set = plot_variables(variables=['E'], xlabel=r'Density (cm$^{-3}$)')
    plot_set_index = 1
    
    for i in range(ion_index, len(atm['Species'])-1):
        t_plot_set.append_variable(atm['Species'][i])
        plot_set_index += 1
        # reach 5 variables for plot or run out of variables (excluding e)
        if (plot_set_index == 5) or (i == len(atm['Species'])-2):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel=r'Density (cm$^{-3}$)')
            plot_set_index = 0
    
    pdf_ion_name = 'density_temp_i.pdf'
    plot_profiles(plot_path + pdf_ion_name, 
                  alt, plot_sets, atm['Profiles'], 
                  plot_title = 'Number Density (Ions)')
    
    # Merge temporary files
    writer = PdfWriter()
    
    reader1 = PdfReader(plot_path + pdf_neutral_name)
    for page in reader1.pages:
        writer.add_page(page)
        
    reader2 = PdfReader(plot_path + pdf_ion_name)
    for page in reader2.pages:
        writer.add_page(page)
    
    with open(plot_path + 'density.pdf', "wb") as fp:
        writer.write(fp)
    
    print(f'Merged {pdf_neutral_name} and {pdf_ion_name} into density.pdf.')

    # Clean up temp files
    os.remove(plot_path + pdf_neutral_name)
    print(f'Removed {pdf_neutral_name}.')
    os.remove(plot_path + pdf_ion_name)
    print(f'Removed {pdf_ion_name}.')
    
#%% Plot mixing ratios


elif p == 2: # Plot mixing ratios
    # Read in atm1D.out
    atm = sgr.atm('../runs/' + run_name + '/output/atm1D.out')
    alt = 1.E-5 * atm['Profiles']['Altitude'] # convert cm to km
    
    # Calculate mixing ratios
    vmr = atm['Profiles'].loc[:,'N_total':'E'].copy()
    for i_mol in range(1, atm['n_mol']+1):
        vmr[atm['Species'][i_mol]] = (vmr[atm['Species'][i_mol]] 
                                      / vmr['N_total'])
    vmr = vmr.drop(columns = 'N_total')
    
    # Plot neutrals and save to temporary file
    plot_sets = []
    t_plot_set = plot_variables(variables=[], xlabel='Mixing Ratio')
    plot_set_index = 0
    for i in range(len(vmr.columns)):
        if 'P' in vmr.columns[i]: # stop when encounter ion
            if t_plot_set.variables: plot_sets.append(t_plot_set)
            ion_index = i
            break
        t_plot_set.append_variable(vmr.columns[i])
        plot_set_index += 1
        if (plot_set_index == 5) or (i == len(vmr.columns)-1):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel='Mixing Ratio')
            plot_set_index = 0
    
    pdf_neutral_name = 'vmr_temp_n.pdf'
    plot_profiles(plot_path + pdf_neutral_name, 
                  alt, plot_sets, vmr,
                  plot_title = 'Mixing Ratio (Neutrals)')
    
    # Plot electrons and ions and save to temporary file
    # move electrons to front
    plot_sets = []
    t_plot_set = plot_variables(variables=['E'], xlabel='Mixing Ratio')
    plot_set_index = 1
    
    for i in range(ion_index, len(vmr.columns)-1):
        t_plot_set.append_variable(vmr.columns[i])
        plot_set_index += 1
        # reach 5 variables for plot or run out of variables (excluding e)
        if (plot_set_index == 5) or (i == len(atm['Species'])-2):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel='Mixing Ratio')
            plot_set_index = 0
            
    pdf_ion_name = 'vmr_temp_i.pdf'
    plot_profiles(plot_path + pdf_ion_name, 
                  alt, plot_sets, vmr,
                  plot_title = 'Mixing Ratio (Ions)')
        
    # Merge temporary files
    writer = PdfWriter()
    
    reader1 = PdfReader(plot_path + pdf_neutral_name)
    for page in reader1.pages:
        writer.add_page(page)
        
    reader2 = PdfReader(plot_path + pdf_ion_name)
    for page in reader2.pages:
        writer.add_page(page)
    
    with open(plot_path + 'vmr.pdf', "wb") as fp:
        writer.write(fp)
    
    print(f'Merged {pdf_neutral_name} and {pdf_ion_name} into vmr.pdf.')


elif p == 3: # Plot single species
    smolec = input('Enter species > ').strip()
    
    # check if species exists
    nmolec = sgr.molecules('../runs/'+run_name+'/input/nmolecules.dat') 
    imolec = sgr.molecules('../runs/'+run_name+'/input/imolecules.dat') 
    molec = pd.concat([nmolec,imolec], ignore_index=True)
    if (molec['Species'] == smolec).any():
        print('Finding and plotting reactions for ' + smolec + '...')
    else: 
        print('Species ' + smolec + ' not found. Exiting...')
        sys.exit()
    
    
    # Plot summary page
    species_summary = sgr.summary('../runs/' + run_name + '/output/molecules/'
                                  + smolec + '.OUT')
    alt = species_summary['Altitude'] # convert cm to km
    
    plot_sets = [plot_variables(variables=['Density'], 
                                xlabel=r'Density (cm$^{-3}$)'),
                 plot_variables(variables=['Mixing Ratio'], 
                                xlabel='Mixing Ratio'),
                 plot_variables(variables=['Production (Net)', 'Loss (Net)', 
                                           'del(flux)' ,'-del(flux)'], 
                                xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)'),
                 plot_variables(variables=['Flux (up)', 'Flux (down)'], 
                                xlabel=r'Flux (cm$^{-2}$ s$^{-1}$)')]
    
    data = species_summary.loc[:,['Density', 'Mixing Ratio', 
                                  'Production (Net)', 'Loss (Net)']]
    
    source, sink = split_up_down(species_summary.loc[:,'cvg_flx'])
    data['del(flux)'] = source
    data['-del(flux)'] = sink
    flux_up, flux_down = split_up_down(species_summary.loc[:,'Flux'])
    flux_up[240] = np.nan
    data['Flux (up)'] = flux_up
    flux_down[240] = np.nan
    data['Flux (down)'] = flux_down

    pdf_sum_name = smolec + '_temp_sum.pdf'
    plot_profiles(plot_path + pdf_sum_name, 
                  alt, plot_sets, data, 
                  plot_title = 'Summary for ' + smolec)
    
    # Get consolidated reaction rates sorted by column rates
    alt, rates, reactions_grouped, sort_order = get_sorted_reactions(run_name)
    # Filter relevant reactions
    reactions_prod = []
    rates_prod = []
    reactions_loss = []
    rates_loss = []
    for i_rct in range(0, len(sort_order)):
        reaction_group = reactions_grouped.get_group(sort_order.index[i_rct])
        # if a product
        if ((reaction_group['Product1'].iloc[0] == smolec) or
            (reaction_group['Product2'].iloc[0] == smolec) or 
            (reaction_group['Product3'].iloc[0] == smolec) or 
            (reaction_group['Product4'].iloc[0] == smolec)):
            reactions_prod.append(reaction_group.iloc[0])
            rates_prod.append(np.sum(rates[:,reaction_group.index], axis=1))
            
        # if a reactant
        if ((reaction_group['Reactant1'].iloc[0] == smolec) or
            (reaction_group['Reactant2'].iloc[0] == smolec)):
            reactions_loss.append(reaction_group.iloc[0])
            rates_loss.append(np.sum(rates[:,reaction_group.index], axis=1))
    
    reactions_prod = pd.DataFrame(reactions_prod)
    rates_prod = np.transpose(rates_prod)
    reactions_loss = pd.DataFrame(reactions_loss)
    rates_loss = np.transpose(rates_loss)

    # Plot production reactions
    data = pd.DataFrame()
    plot_sets = []
    t_plot_set = plot_variables(variables=[], 
                                xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
    plot_set_index = 0
    for i_rct in range(len(reactions_prod)):
        title = reaction_title(reactions_prod['Reactant1'].iloc[i_rct],
                               reactions_prod['Reactant2'].iloc[i_rct],
                               reactions_prod['Product1'].iloc[i_rct],
                               reactions_prod['Product2'].iloc[i_rct],
                               reactions_prod['Product3'].iloc[i_rct],
                               reactions_prod['Product4'].iloc[i_rct])
        t_plot_set.append_variable(title)
        data[title] = rates_prod[:,i_rct]
        plot_set_index += 1
        if (plot_set_index == 5) or (i_rct == len(reactions_prod)-1):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
            plot_set_index = 0
    
    pdf_prod_name = smolec + '_temp_prod.pdf'
    plot_profiles(plot_path + pdf_prod_name, 
                  alt, plot_sets, data, 
                  plot_title = 'Production Reactions')
    
    # Plot loss reactions
    data = pd.DataFrame()
    plot_sets = []
    t_plot_set = plot_variables(variables=[], 
                                xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
    plot_set_index = 0
    for i_rct in range(len(reactions_loss)):
        title = reaction_title(reactions_loss['Reactant1'].iloc[i_rct],
                               reactions_loss['Reactant2'].iloc[i_rct],
                               reactions_loss['Product1'].iloc[i_rct],
                               reactions_loss['Product2'].iloc[i_rct],
                               reactions_loss['Product3'].iloc[i_rct],
                               reactions_loss['Product4'].iloc[i_rct])
        t_plot_set.append_variable(title)
        data[title] = rates_loss[:,i_rct]
        plot_set_index += 1
        if (plot_set_index == 5) or (i_rct == len(reactions_loss)-1):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
            plot_set_index = 0
    
    pdf_loss_name = smolec + '_temp_loss.pdf'
    plot_profiles(plot_path + pdf_loss_name, 
                  alt, plot_sets, data, 
                  plot_title = 'Loss Reactions')
    
    # Merge temporary files
    writer = PdfWriter()
    
    reader1 = PdfReader(plot_path + pdf_sum_name)
    for page in reader1.pages:
        writer.add_page(page)
        
    reader2 = PdfReader(plot_path + pdf_prod_name)
    for page in reader2.pages:
        writer.add_page(page)
        
    reader3 = PdfReader(plot_path + pdf_loss_name)
    for page in reader3.pages:
        writer.add_page(page)

    with open(plot_path + smolec + '.pdf', "wb") as fp:
        writer.write(fp)
    
    print(f'Merged {pdf_sum_name}, {pdf_prod_name}, and {pdf_loss_name} ' + 
          f'into {smolec}.pdf.')

    # Clean up temp files
    os.remove(plot_path + pdf_sum_name)
    print(f'Removed {pdf_sum_name}.')
    os.remove(plot_path + pdf_prod_name)
    print(f'Removed {pdf_prod_name}.')
    os.remove(plot_path + pdf_loss_name)
    print(f'Removed {pdf_loss_name}.')
    
    
elif p == 4: # Plot all reactions
    # Get consolidated reaction rates sorted by column rates
    alt, rates, reactions_grouped, sort_order = get_sorted_reactions(run_name)
    reactions_sorted = []
    rates_sorted = []
    for i_rct in range(0, len(sort_order)):
        reaction_group = reactions_grouped.get_group(sort_order.index[i_rct])
        reactions_sorted.append(reaction_group.iloc[0])
        rates_sorted.append(np.sum(rates[:,reaction_group.index], axis=1))
    
    reactions_sorted = pd.DataFrame(reactions_sorted)
    rates_sorted = np.transpose(rates_sorted)

    # Plot reactions
    data = pd.DataFrame()
    plot_sets = []
    t_plot_set = plot_variables(variables=[], 
                                xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
    plot_set_index = 0
    for i_rct in range(len(reactions_sorted)):
        title = reaction_title(reactions_sorted['Reactant1'].iloc[i_rct],
                               reactions_sorted['Reactant2'].iloc[i_rct],
                               reactions_sorted['Product1'].iloc[i_rct],
                               reactions_sorted['Product2'].iloc[i_rct],
                               reactions_sorted['Product3'].iloc[i_rct],
                               reactions_sorted['Product4'].iloc[i_rct])
        t_plot_set.append_variable(title)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", 
                                  category=pd.errors.PerformanceWarning)
            data[title] = rates_sorted[:,i_rct]
        plot_set_index += 1
        if (plot_set_index == 5) or (i_rct == len(reactions_sorted)-1):
            plot_sets.append(t_plot_set)
            t_plot_set = plot_variables(variables=[], 
                                        xlabel=r'Rate (cm$^{-3}$ s$^{-1}$)')
            plot_set_index = 0
    
    pdf_name = 'rates.pdf'
    plot_profiles(plot_path + pdf_name, 
                  alt, plot_sets, data, 
                  plot_title = 'Reaction Rates')