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

def covmat_plots(label, fp1_table, fp2_table):
    # Separating data sets to produce a separate covmat for each one
    datasets = fp1_table.dataset.unique()
    for ds in datasets:
        fp1_table_ds = fp1_table.query(f'dataset == "{ds}"')
        fp2_table_ds = fp2_table.query(f'dataset == "{ds}"')

        T_fp2 = fp2_table_ds["theory_central"]
        D = fp1_table_ds["data_central"]

        # Calculating errors on T_fp1 by taking standard deviation
        T_fp1_reps = fp1_table_ds.loc[:, fp1_table.columns.str.contains("rep")]
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
        cbar.set_label(label="% of data", fontsize=12)
        cbar.ax.tick_params(labelsize=12)
        ax.set_title(f"{ds}", fontsize=15)
        plt.savefig(f"../../plots/covmats/covmats_{label}_{ds}.png")

        # Diag element plot
        sqrtdiags = np.sqrt(np.diag(normcovmat))
        fig2, ax2 = plt.subplots(figsize=(15,6))
        ax2.plot(100*sqrtdiags, '-o', color="darkorange")
        ax2.set_ylabel("% of data")
        ax2.set_title(f"Diagonal elements of covariance matrix for {ds}",
                      fontsize=15)
        plt.savefig(f"../../plots/covmats/diag_covmat_{label}_{ds}.png")
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
