# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 01:00:14 2024

@author: hsrivast
"""

import os

# Specify the path to the new working directory
new_working_directory='X:/fast/gilbert_p/fg_data/SEATRAC/TB_hackday_2023/seatrac-hackday-2023/liu_etal/Eigenplot'

# Change the current working directory to the new directory
os.chdir(new_working_directory)

import matplotlib as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import itertools
import os
from os.path import join as opj
import sys
from pngpdf import PngPdfPages
from biplot import biplot
from myboxplot import swarmbox

sns.set_style('whitegrid')

import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
GSE218157_family_soft = pd.read_csv("X:/fast/gilbert_p/fg_data/SEATRAC/TB_hackday_2023/seatrac-hackday-2023/liu_etal/GSE218157_family.soft.csv")

# Convert the DataFrame to a pandas DataFrame (if not already in pandas format)
meta_data = GSE218157_family_soft

# Set the row names of the DataFrame to the 'description' column
meta_data.index = meta_data['description']

# Recode the 'time_after_bcg' column into a new 'day' column
meta_data['day'] = meta_data['time_after_bcg'].map({
    "pre": "0", 
    "d2": "2", 
    "wk2": "14", 
    "wk12": "84"
})





eg_fn = "module_eigen_genes.csv"
eg_ = pd.read_csv(eg_fn).rename({'Unnamed: 0':'description'}, axis=1)
eg = eg_.drop(eg_.columns[-1], axis=1)
eg.set_index(eg.columns[0], inplace=True)


merged_data = pd.merge(eg, meta_data[['vax_group',"day","pitt_id"]], left_index=True, right_index=True)

module_names = eg.columns.tolist()
modules = module_names[1:]

treatments = meta_data['vax_group'].unique()


with PngPdfPages('X:/fast/gilbert_p/fg_data/SEATRAC/TB_hackday_2023/seatrac-hackday-2023/liu_etal/Eigenplot/eigen_gene_plots.pdf') as pdf:
    for trt in treatments:
        for m in modules:
            figh = plt.figure(figsize=(15,10))
            plotdf = merged_data.loc[merged_data['vax_group'] == trt]
            swarmbox(x='day', y=m, connect=True, connect_on=['pitt_id'], data=plotdf, order=['0', '2', '14', '84'])
            plt.title(f'{trt}')
            pdf.savefig(figh)


with PngPdfPages('X:/fast/gilbert_p/fg_data/SEATRAC/TB_hackday_2023/seatrac-hackday-2023/liu_etal/Eigenplot/eigen_gene_plots_pooled.pdf') as pdf:
    for m in modules:
        figh = plt.figure(figsize=(15,10))
        sns.lineplot(x='day', y=m, hue='vax_group', data=merged_data)
        plt.legend(loc='upper left', bbox_to_anchor=(1,1))
        pdf.savefig(figh)

