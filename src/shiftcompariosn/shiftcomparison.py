#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:15:00 2020

@author: s1303034
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm, colors as mcolors
import sys

np.set_printoptions(threshold=sys.maxsize)

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
    # Renaming ticklabels
    ticklabel_dict = {"NMCPD": "NMC",
                      "NMCPD_D": "NMC",
                      "SLACD": "SLAC",
                      "BCDMSD": "BCDMS",
                      "DYE886_D": "NuSea"}
    ticklabels = [ticklabel_dict[ticklabel] for ticklabel in ticklabels]
    startlocs = [unique_ds[x][1] for x in range(len(unique_ds))]
    startlocs += [len(labels)]
    ticklocs = [0 for x in range(len(startlocs) - 1)]
    for i in range(len(startlocs) - 1):
        ticklocs[i] = 0.5 * (startlocs[i + 1] + startlocs[i])
    return ticklocs, ticklabels, startlocs


fp1_table = pd.read_table(
    "../observables/fp1_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

fp2_table = pd.read_table(
    "../observables/fp2_global/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)



T_fp1 = fp1_table["theory_central"]
T_fp2 = fp2_table["theory_central"]
shift = T_fp2.values - T_fp1.values

# Calculating errors on T_fp1 by taking standard deviation
T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]
T_fp1_errs = T_fp1_reps.std(axis=1)

fig, ax = plt.subplots(figsize=(20, 7))
ax.errorbar(range(len(T_fp1)),
                (shift / T_fp1).values,
                yerr= (T_fp1_errs / T_fp1).values,
                marker=".",
                ls = '', label=r"$(T_i^d[f_d] - T_i^d[f_p])$")
ticklocs, ticklabels, startlocs = matrix_plot_labels(fp1_table)
plt.xticks(ticklocs, ticklabels, rotation=40, fontsize=16)
# Shift startlocs elements 0.5 to left so lines are between indexes
startlocs_lines = [x - 0.5 for x in startlocs]
ax.margins(x=0, y=0)
ax.set_title(r"Normalised to $T_i^d[f_p]$", fontsize=28)
ax.set_ylim([-0.25,0.25])
plt.yticks(fontsize=16)
ymin, ymax = ax.get_ylim()
xmin, xmax = ax.get_xlim()
ax.hlines(0, xmin, xmax, linestyles="-")
ax.vlines(startlocs_lines, ymin, ymax, linestyles="dashed")
ax.legend(fontsize="20")
    
    


