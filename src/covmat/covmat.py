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
    explabels = [x[0] for x in df.index]
    datasetlabels = [x[1] for x in df.index]
    points = [x[2] for x in df.index]
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

def matrix_plot(matrix, labeldf, descrip, label):
    fig, ax = plt.subplots(figsize=(6,6))
    matrixplot = ax.matshow(100*matrix,
                            cmap=cm.Spectral_r)
    cbar=fig.colorbar(matrixplot, fraction=0.046, pad=0.04)
    cbar.set_label(label="% of data", fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    ax.set_title(f"Covariance matrix: {descrip} {label}", fontsize=15)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(labeldf)
    plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=15)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(ticklocs, ticklabels, fontsize=15)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.vlines(startlocs_lines, 0, 1, transform=ax.get_xaxis_transform(), linestyles="dashed")
    ax.hlines(startlocs_lines,0, 1, transform=ax.get_yaxis_transform(), linestyles='dashed')
    ax.margins(x=0, y=0)
    plt.savefig(f"../../plots/covmats/covmats_{descrip}_{label}.png")
    return fig

def covmat_plots(label, fp1_table, fp2_table, fp1_covmat):
    # Separating data sets to produce a separate covmat for each one
#    datasets = fp1_table.dataset.unique()

   # # Cutting down experimental covmat
   # expcovmat = ((fp1_covmat.values[3:]).T[3:].T).astype(np.float64)
    expcovmat = fp1_covmat.values

    T_fp2 = fp2_table["theory_central"]
    D = fp1_table["data_central"]

    # Calculating errors on T_fp1 by taking standard deviation
    T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]
    nrep = len(T_fp1_reps.values.T)

    T_fp2_repeat = np.tile(T_fp2.values, (nrep,1)).T
    deltas = T_fp1_reps.values - T_fp2_repeat

    covmat = (1/nrep) * deltas@deltas.T
    normcovmat = covmat/np.outer(D.values, D.values)
    normexpcovmat = expcovmat/np.outer(D.values, D.values)
    expsqrtdiags = np.sqrt(np.diag(normexpcovmat))
    impactcovmat = (covmat + expcovmat)/expcovmat

    fig_th = matrix_plot(normcovmat, fp1_table, "theory", label)
    fig_imp = matrix_plot(impactcovmat, fp1_table, "impact", label)

    # Diag element plot
    sqrtdiags = np.sqrt(np.diag(normcovmat))
    fig_diag, ax2 = plt.subplots(figsize=(15,6))
    ax2.plot(100*sqrtdiags, '-o', color="darkorange", label="S")
    ax2.plot(100*expsqrtdiags, '-o', color="purple", label="C")
    ax2.set_ylabel("% of data", fontsize=15)
    ax2.set_title(f"Diagonal elements of covariance matrix: {label}",
                  fontsize=15)
    ticklocs, ticklabels, startlocs = matrix_plot_labels(fp1_table)
    plt.xticks(ticklocs, ticklabels, rotation=30, ha="right", fontsize=15)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ymin, ymax = ax2.get_ylim()
    xmin, xmax = ax2.get_xlim()
    ax2.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
    ax2.margins(x=0, y=0)
    ax2.legend(fontsize=15)
    plt.savefig(f"../../plots/covmats/diag_covmat_{label}.png")

    return fig_th, fig_imp, fig_diag

# Loading DIS and global experiment tables

fp1_table_DIS = pd.read_table(
    "../observables/fp1/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

fp2_table_DIS = pd.read_table(
    "../observables/fp2/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

fp1_table_global = pd.read_table(
    "../observables/fp1_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

fp2_table_global = pd.read_table(
    "../observables/fp2_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

# Loading DIS and experiment covmats

fp1_covmat_DIS = pd.read_table(
    "../observables/fp1/output/tables/experiments_covmat.csv",
    dtype={"user_id": float},
    index_col=[0,1,2], header=[0,1,2]
)

fp2_covmat_DIS = pd.read_table(
    "../observables/fp2/output/tables/experiments_covmat.csv",
    dtype={"user_id": float},
    index_col=[0,1,2], header=[0,1,2]
)

fp1_covmat_global = pd.read_table(
    "../observables/fp1_global/output/tables/experiments_covmat.csv",
    dtype={"user_id": float},
    index_col=[0,1,2], header=[0,1,2]
)

fp2_covmat_global = pd.read_table(
    "../observables/fp2_global/output/tables/experiments_covmat.csv",
    dtype={"user_id": float},
    index_col=[0,1,2], header=[0,1,2]
)

# Plotting

covmat_plots("DIS", fp1_table_DIS, fp2_table_DIS,
             fp1_covmat_DIS)

covmat_plots("global_proton", fp1_table_global, fp2_table_global,
             fp1_covmat_global)
