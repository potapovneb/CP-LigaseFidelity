#!/usr/bin/env python3

###############################################################################
# Ligase Fidelity Profiling
# Copyright (C) 2022 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('ligation_data_dir')
args = parser.parse_args()


def plot_matrix(input_file):
    # Load data
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
                        cbar=False,
                        cbar_kws={"shrink": .5},
                        annot=False,
                        xticklabels=True,
                        yticklabels=True,
                        ax=ax)

    plt.xlabel('')
    plt.ylabel('')

    # Save the figure
    output_file = input_file.replace('.csv','.png')
    plt.savefig(output_file)


def plot_fidelity(input_file):
    # Load data
    fidelity = pd.read_csv(input_file, sep=',')

    # Normalize data to 100,000
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
    # ax.set(xlim=(0, 800), ylabel="", xlabel="Ligation events")
    ax.set(ylabel="", xlabel="Ligation events")
    sns.despine(left=True, bottom=True)
    
    # Save the figure
    output_file = input_file.replace('.csv','.png')
    plt.savefig(output_file)


def plot_mismatch(input_file):
    # Load data
    mismatch = pd.read_csv(input_file, sep=',')
    mismatch['Mismatch'] = mismatch.apply(lambda row: row['Bottom'] + row['Top'],axis=1)
    mismatch = mismatch.sort_values(by='Count',ascending=False)

    # Calculate percentages
    total = mismatch['Count'].sum()
    mismatch['Count'] = mismatch['Count'] / total * 100

    # Discard Watson-Crick pairs from the plot
    mismatch = mismatch[(mismatch['Mismatch'] != 'AT') & (mismatch['Mismatch'] != 'TA') & (mismatch['Mismatch'] != 'CG') & (mismatch['Mismatch'] != 'GC')]
    mismatch = mismatch[['Mismatch','Count']]

    # Set figure style
    # sns.set_theme(style="whitegrid",font="courier new")
    sns.set_theme(style="whitegrid")

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(7, 7))

    # Define colors
    col = sns.color_palette("Set2")

    # Plot the total number of ligation events
    sns.barplot(x='Mismatch', y='Count', data=mismatch, color=col[1])
    # ax.set(ylim=(0,5), ylabel='Mismatch frequency, %')
    ax.set(ylabel='Mismatch frequency, %')
    sns.despine(left=True, bottom=True)

    # Save the figure
    output_file = input_file.replace('.csv','.png')
    plt.savefig(output_file)

plot_matrix(args.ligation_data_dir   + '/' + '06_matrix.csv')
plot_fidelity(args.ligation_data_dir + '/' + '07_fidelity.csv')
plot_mismatch(args.ligation_data_dir + '/' + '08_mismatch-e.csv')
plot_mismatch(args.ligation_data_dir + '/' + '09_mismatch-m.csv')
