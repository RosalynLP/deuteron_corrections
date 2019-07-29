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
    """Returns the tick locations and labels, and the starting
    point values for each category,  based on a dataframe
    to be plotted. The dataframe is assumed to be multiindexed by
    (process, dataset, points) or else (dataset, points). The tick
    location is in the centre of the dataset, and labelling is by
    the outermost index of the multiindex."""
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


fp1_table = pd.read_table(
    "./fp1/output/tables/experiment_result_table.csv", dtype={"user_id": float}
)

fp2_table = pd.read_table(
    "./fp2/output/tables/experiment_result_table.csv", dtype={"user_id": float}
)

# Cutting to same size for comparison
#fp1_table = fp1_table.query('dataset == "BCDMSD"')
#fp2_table = fp2_table.query('dataset == "BCDMSD"')

T_fp1 = fp1_table["theory_central"]
T_fp2 = fp2_table["theory_central"]

# Calculating errors on T_fp1 by taking standard deviation
T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]
T_fp1_errs = T_fp1_reps.std(axis=1)

# Plotting
fig, ax = plt.subplots(figsize=(15, 6))
ax.errorbar(range(len(T_fp1)),
            T_fp1 / T_fp2,
            yerr=T_fp1_errs / T_fp2,
            marker=".",
            ls = '')
ticklocs, ticklabels, startlocs = matrix_plot_labels(fp1_table)
plt.xticks(ticklocs, ticklabels, rotation=45, fontsize=15)
# Shift startlocs elements 0.5 to left so lines are between indexes
startlocs_lines = [x - 0.5 for x in startlocs]
ymin, ymax = ax.get_ylim()
xmin, xmax = ax.get_xlim()
ax.hlines(1, xmin, xmax, linestyles="-")
ax.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
ax.margins(x=0, y=0)
ax.set_title(r"$T_i^d[f_d] / \langle T_i^d[f_p] \rangle$", fontsize=20)
plt.savefig("../../plots/observables/observable_ratio.png")
