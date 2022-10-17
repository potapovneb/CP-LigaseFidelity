#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('ligation_data_dir')
args = parser.parse_args()

def plot_matrix(input_file):
    data = pd.read_csv(input_file, sep=',', index_col=0)

    # normalize ligation events per 100,0000
    data_log = data / data.to_numpy().sum() * 1e5

    # replace 0 with NANs (so log10() works fine)
    data_log = data_log.replace(0, np.nan)
    data_log = np.log10(data_log)

    # define fixed font for overhang labels
    # sns.set_theme(font="courier new")

    # define colormap
    cmap = sns.color_palette("light:b", as_cmap=True)

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(50, 50))

    # Draw the heatmap
    heatmap = sns.heatmap(data_log,
                        cmap=cmap,
                        square=True,
                        linecolor='white',
                        linewidths=0.1,
                        cbar=True,
                        cbar_kws={"shrink": .5},
                        annot=False,
                        xticklabels=True,
                        yticklabels=True,
                        ax=ax)

    plt.xlabel('')
    plt.ylabel('')

    output_file = input_file.split('.')[0] + '.png'
    plt.savefig(output_file)

def plot_ligation_fidelity(input_file):
    # Load data
    fidelity = pd.read_csv(input_file, sep=',')

    # Normalize data
    total = fidelity['Total'].sum()
    fidelity['Total'] = fidelity['Total'] / total * 1e5
    fidelity['Mismatch'] = fidelity['Mismatch'] / total * 1e5
    fidelity['Correct'] = fidelity['Correct'] / total * 1e5

    # Set figure style
    # sns.set_theme(style="whitegrid",font="courier new")
    sns.set_theme(style="whitegrid")

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(10, 50))

    # Define colors
    col = sns.color_palette("Set2")

    # Plot the total number of ligation events
    sns.barplot(x='Total', y='Overhang', data=fidelity, label='Mismatch', color=col[1])

    # Plot the correct number of ligation events
    sns.barplot(x='Correct', y='Overhang', data=fidelity, label='Correct', color=col[2])

    # Add a legend and informative axis label
    ax.legend(ncol=2, loc="lower right", frameon=True)
    ax.set(xlim=(0, 800), ylabel="", xlabel="Ligation events")
    sns.despine(left=True, bottom=True)
    
    output_file = input_file.split('.')[0] + '.png'
    plt.savefig(output_file)

def plot_mismatch(input_file):
    # Load data
    mismatch = pd.read_csv(input_file, sep=',')
    mismatch['Mismatch'] = mismatch.apply(lambda row: row['Bottom'] + row['Top'],axis=1)
    mismatch = mismatch.sort_values(by='Count',ascending=False)

    total = mismatch['Count'].sum()
    mismatch['Count'] = mismatch['Count'] / total * 100

    mismatch = mismatch[(mismatch['Mismatch'] != 'AT') & (mismatch['Mismatch'] != 'TA') & (mismatch['Mismatch'] != 'CG') & (mismatch['Mismatch'] != 'GC')]
    mismatch = mismatch[['Mismatch','Count']]

    # Set figure style
    # sns.set_theme(style="whitegrid",font="courier new")
    sns.set_theme(style="whitegrid")

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(5, 5))

    # Define colors
    col = sns.color_palette("Set2")

    # Plot the total number of ligation events
    sns.barplot(x='Mismatch', y='Count', data=mismatch, color=col[1])
    ax.set(ylim=(0,5), ylabel='Mismatch frequency')
    sns.despine(left=True, bottom=True)
    
    output_file = input_file.split('.')[0] + '.png'
    plt.savefig(output_file)

plot_matrix(args.ligation_data_dir + '/06_matrix.csv')
plot_ligation_fidelity(args.ligation_data_dir + '/07_fidelity.csv')
plot_mismatch(args.ligation_data_dir + '/08_mismatch-e.csv')
plot_mismatch(args.ligation_data_dir + '/09_mismatch-m.csv')
