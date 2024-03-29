#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:04:39 2019

@author: rosalyn
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    print(ticklabels)
    # Renaming ticklabels
    ticklabel_dict = {"SLACD": "SLAC",
                      "SLACD_dw_ite": "SLAC",
                      "BCDMSD": "BCDMS",
                      "BCDMSD_dw_ite": "BCDMS",
                      "NMCPD": "NMC",
                      "NMCPD_D": "NMC",
                      "NMCPD_dw_ite": "NMC",
                      "DYE886_D": "NuSea",
                      "DYE906_D": "SeaQuest"}
    ticklabels = [ticklabel_dict[ticklabel] for ticklabel in ticklabels]
    startlocs = [unique_ds[x][1] for x in range(len(unique_ds))]
    startlocs += [len(labels)]
    ticklocs = [0 for x in range(len(startlocs) - 1)]
    for i in range(len(startlocs) - 1):
        ticklocs[i] = 0.5 * (startlocs[i + 1] + startlocs[i])
    return ticklocs, ticklabels, startlocs

def plot_observable_ratio(label, fp1_table, fp2_table):

    T_fp1 = fp1_table["theory_central"]
    T_fp2 = fp2_table["theory_central"]

    # Calculating errors on T_fp1 by taking standard deviation
    T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]
    T_fp1_errs = T_fp1_reps.std(axis=1)

    # Plotting
    fig, ax = plt.subplots(figsize=(20, 7))
    ax.errorbar(range(len(T_fp1)),
                T_fp1 / T_fp2,
                yerr=T_fp1_errs / T_fp2,
                marker=".",
                ls = '')
    ticklocs, ticklabels, startlocs = matrix_plot_labels(fp1_table)
    plt.xticks(ticklocs, ticklabels, rotation=40, fontsize=16)
    # Shift startlocs elements 0.5 to left so lines are between indexes
    startlocs_lines = [x - 0.5 for x in startlocs]
    ax.margins(x=0, y=0)
    ax.set_title(r"$T_i^d[f_d^0] / \langle T_i^d[f_s^0] \rangle$", fontsize=28)
    ax.set_ylim([0.85,1.15])
    plt.yticks(fontsize=16)
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.hlines(1, xmin, xmax, linestyles="-")
    ax.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
    plt.savefig(f"../../plots/observables/observable_ratio_{label}.png")
    return fig


# Loading DIS and global experiment tables

fp1_table_global = pd.read_table(
    "./fp1_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

fp2_table_global = pd.read_table(
    "./fp2_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

# Loading iteration fit

fp2_table_global_iteration1 = pd.read_table(
    "./fp2_global_iteration1/output/tables/experiment_result_table.csv",
    dtype={"user_id": float}
)

fp2_table_global_iteration2 = pd.read_table(
    "./fp2_global_iteration2/output/tables/group_dataset_inputs_by_experiment0_group_result_table.csv",
    dtype={"user_id": float}
)

fp1_table_global_iteration2 = pd.read_table(
    "./fp1_global_iteration2/output/tables/group_dataset_inputs_by_experiment0_group_result_table.csv",
    dtype={"user_id": float}
)

# Plotting

plot_observable_ratio("global_proton", fp1_table_global, fp2_table_global)
plot_observable_ratio("ite_1", fp1_table_global, fp2_table_global_iteration1)
plot_observable_ratio("ite_2", fp1_table_global_iteration2, fp2_table_global_iteration2)


