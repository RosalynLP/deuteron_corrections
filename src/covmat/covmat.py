#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:11:05 2019

@author: rosalyn
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

sns.set()

def matrix_plot_labels(df):
    datasetlabels = [x for x in df["dataset"]]
    points = [x for x in df["id"]]
    labels = datasetlabels
    unique_ds = []
    unique_ds.append([labels[0], points[0]])
    for x in range(len(labels) - 1):
        if labels[x + 1] != labels[x]:
            unique_ds.append([labels[x + 1], x + 1])
    ticklabels = [unique_ds[x][0] for x in range(len(unique_ds))]
    startlocs = [unique_ds[x][1] for x in range(len(unique_ds))]
    startlocs += [len(labels)]
    ticklocs = [0 for x in range(len(startlocs) - 1)]
    for i in range(len(startlocs) - 1):
        ticklocs[i] = 0.5 * (startlocs[i + 1] + startlocs[i])
    return ticklocs, ticklabels, startlocs

def covmat_plots(label, fp1_table, fp2_table):
    # Separating data sets to produce a separate covmat for each one
#    datasets = fp1_table.dataset.unique()
    
    T_fp2 = fp2_table["theory_central"]
    D = fp1_table["data_central"]

    # Calculating errors on T_fp1 by taking standard deviation
    T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]
    nrep = len(T_fp1_reps.values.T)

    T_fp2_repeat = np.tile(T_fp2.values, (nrep,1)).T
    deltas = T_fp1_reps.values - T_fp2_repeat

    covmat = (1/nrep) * deltas@deltas.T
    normcovmat = covmat/np.outer(D.values, D.values)

    # Full covmat plot
    fig, ax = plt.subplots(figsize=(6,6))
    matrixplot = ax.matshow(100*normcovmat,
                            cmap=cm.Spectral_r)
    cbar=fig.colorbar(matrixplot, fraction=0.046, pad=0.04)
    cbar.set_label(label="% of data", fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    ax.set_title(f"Covariance matrix: {label}", fontsize=15)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(fp1_table)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=15)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels, fontsize=15)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
    ax.hlines(startlocs_lines, xmin, xmax, linestyles='dashed')
    ax.margins(x=0, y=0)
    plt.savefig(f"../../plots/covmats/covmats_{label}.png")

    # Diag element plot
    sqrtdiags = np.sqrt(np.diag(normcovmat))
    fig2, ax2 = plt.subplots(figsize=(15,6))
    ax2.plot(100*sqrtdiags, '-o', color="darkorange")
    ax2.set_ylabel("% of data", fontsize=15)
    ax2.set_title(f"Diagonal elements of covariance matrix: {label}",
                  fontsize=15)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right", fontsize=15)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ymin, ymax = ax2.get_ylim()
    xmin, xmax = ax2.get_xlim()
    ax2.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
    ax2.margins(x=0, y=0)
    plt.savefig(f"../../plots/covmats/diag_covmat_{label}.png")
    
    return fig, fig2

# Loading DIS and global experiment tables

fp1_table_DIS = pd.read_table(
    "../observables/fp1/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

fp2_table_DIS = pd.read_table(
    "../observables/fp2/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

#fp1_table_global = pd.read_table(
#    "../observables/fp1_global/output/tables/experiment_result_table.csv",
#    dtype={"user_id": float}
#)

fp2_table_global = pd.read_table(
    "../observables/fp2_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

# Plotting

covmat_plots("DIS", fp1_table_DIS, fp2_table_DIS)

covmat_plots("global_proton", fp1_table_DIS, fp2_table_global)
