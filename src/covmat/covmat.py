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

fp1_table = pd.read_table(
    "../observables/fp1/output/tables/experiment_result_table.csv", dtype={"user_id": float}
)

fp2_table = pd.read_table(
    "../observables/fp2/output/tables/experiment_result_table.csv", dtype={"user_id": float}
)

# Cutting to same size for comparison
#fp1_table = fp1_table.query('dataset == "BCDMSD"')
#fp2_table = fp2_table.query('dataset == "BCDMSD"')

T_fp1 = fp1_table["theory_central"]
T_fp2 = fp2_table["theory_central"]

# Calculating errors on T_fp1 by taking standard deviation
T_fp1_reps = fp1_table.loc[:, fp1_table.columns.str.contains("rep")]

nrep = len(T_fp1_reps.values.T)

T_fp2_repeat = np.tile(T_fp2.values, (nrep,1)).T
deltas = T_fp1_reps.values - T_fp2_repeat

covmat = (1/nrep) * deltas@deltas.T

fig, ax = plt.subplots(figsize=(6,6))
ax.matshow(covmat, cmap=cm.Spectral_r)