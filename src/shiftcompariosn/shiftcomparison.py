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
import sys

np.set_printoptions(threshold=sys.maxsize)

sns.set_style(style="white")


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

cfacnum_table = pd.read_table(
    "../cfactor/deuteron1/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

cfacdenom_table = pd.read_table(
    "../cfactor/deuteron0/output/tables/experiment_result_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

fp1_kintable = pd.read_table(
    "../observables/fp1_global/output/tables/experiments_xq2_table.csv",
    dtype={"user_id": float},
    index_col=[0,1,2]
)

T_fp1 = fp1_table["theory_central"]
T_fp2 = fp2_table["theory_central"]
shift = T_fp2.values - T_fp1.values
cfac = T_fp1.values/T_fp2.values
    

# MMHT model - generate cfac
def generate_MMHT_cfac(xvals, N, c1, c2, c3e8, xp):
    c1 = np.full_like(xvals, c1)
    c2 = np.full_like(xvals, c2)
    c3 = np.full_like(xvals, c3e8/1e8)
    xp = np.full_like(xvals, xp)
 #   c_3params = (1 + 0.01*N)*(1 + 0.01*c1*(np.log(xp/xvals))**2)
    c_4params = (1 + 0.01*N)*(1 + 0.01*c2*(np.log(xvals/xp))**2
                              + 0.01*c3*(np.log(xvals/xp))**20)
    return c_4params

c_4params = generate_MMHT_cfac(fp1_kintable["x"], 0.589, -0.116, -0.384, 0.0489, 0.03)
#c_3params = generate_MMHT_cfac(fp1_kintable["x"],-0.490, 0.349, -0.444, 3.40, 0.05)

fig, ax = plt.subplots(figsize=(8,6))

plt.plot(fp1_kintable["x"], cfac, "o", label="from fitted deuteron PDF")
#plt.plot(fp1_kintable["x"], c_3params, "o", label="MMSTWW 3 params")
plt.plot(fp1_kintable["x"], c_4params, "o", label="MMHT2014 NNLO 4 params")
ax.set_xlim(0.01,1)
ax.set_xscale("log")
#ax.set_ylim([0.94,1.04])
xmin, xmax = ax.get_xlim()
ax.hlines(1, xmin, xmax, linestyles="-")
ax.legend()


