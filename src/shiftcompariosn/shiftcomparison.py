#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:15:00 2020

@author: s1303034
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

np.set_printoptions(threshold=sys.maxsize)

sns.set_style(style="white")

import apfel

xlha = np.logspace(-2, 0, num=30)[:-2]

# activate some options
apfel.SetPerturbativeOrder(2)
apfel.SetMassScheme("FONLL-C")

# Setting proton only PDF set
apfel.SetPDFSet("NNPDF31_nnlo_as_0118_proton_only")

# initializes integrals on the grids
apfel.InitializeAPFEL_DIS()

eps = 1e-10

Q0 = 2**(1/2) - eps
Q =  10**(1/2)

apfel.ComputeStructureFunctionsAPFEL(Q0,Q)

F2s_proton = []
for x in xlha:
    F2s_proton.append(apfel.F2total(x))
    
# Setting deuteron only PDF
apfel.SetPDFSet("NNPDF31_nnlo_as_0118_deuteron_only")
apfel.InitializeAPFEL_DIS()
apfel.ComputeStructureFunctionsAPFEL(Q0,Q)

F2s_deuteron = []
for x in xlha:
    F2s_deuteron.append(apfel.F2total(x))
    
# Nuclear fit deuteron PDF

apfel.SetPDFSet("nNNPDF20_nlo_as_0118_D2")
apfel.InitializeAPFEL_DIS()
apfel.ComputeStructureFunctionsAPFEL(Q0,Q)

F2s_nuc = []
for x in xlha:
    F2s_nuc.append(apfel.F2total(x))
    
cfac_deut = [a/b for (a, b) in zip(F2s_deuteron, F2s_proton)]
cfac_nuc = [a/b for (a, b) in zip(F2s_nuc, F2s_proton)]
      
    
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

c_4params = generate_MMHT_cfac(xlha, 0.589, -0.116, -0.384, 0.0489, 0.03)
#c_3params = generate_MMHT_cfac(fp1_kintable["x"],-0.490, 0.349, -0.444, 3.40, 0.05)


fig, ax = plt.subplots(figsize=(8,6))

plt.plot(xlha, cfac_deut, label="from fitted deuteron PDF")
plt.plot(xlha, cfac_nuc, label="from nNNPDF2.0")
#plt.plot(fp1_kintable["x"], c_3params, "o", label="MMSTWW 3 params")
plt.plot(xlha, c_4params, "-", label="MMHT2014 NNLO 4 params")
ax.set_xlim(0.01,1)
ax.set_xscale("log")
#ax.set_ylim([0.94,1.04])
xmin, xmax = ax.get_xlim()
ax.set_xlabel("x")
ax.set_ylabel("Correction factor")
plt.title(r"$Q^2$ = 10 GeV$^2$")
ax.hlines(1, xmin, xmax, linestyles="-")
ax.legend()


