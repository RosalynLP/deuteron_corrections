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

sns.set_style(style="ticks")

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


F2s_proton = []
for x in xlha:
    F2reps = []
    for i in range(1,101):
        apfel.SetReplica(i)
        apfel.ComputeStructureFunctionsAPFEL(Q0,Q)
        F2reps.append(apfel.F2total(x))
    F2s_proton.append(F2reps)
    
F2s_proton = np.array(F2s_proton)    
errs_proton = np.std(F2s_proton, axis=1)
cvs_proton = np.mean(F2s_proton, axis=1)

# Setting deuteron only PDF
apfel.SetPDFSet("NNPDF31_nnlo_as_0118_deuteron_only")
apfel.InitializeAPFEL_DIS()

F2s_deuteron = []
for x in xlha:
    F2reps = []
    for i in range(1,101):
        apfel.SetReplica(i)
        apfel.ComputeStructureFunctionsAPFEL(Q0,Q)
        F2reps.append(apfel.F2total(x))
    F2s_deuteron.append(F2reps)
        
F2s_deuteron = np.array(F2s_deuteron)    
errs_deuteron = np.std(F2s_deuteron, axis=1)
cvs_deuteron = np.mean(F2s_deuteron, axis=1)

# Nuclear fit deuteron PDF

apfel.SetPDFSet("nNNPDF20_nlo_as_0118_D2")
apfel.InitializeAPFEL_DIS()

F2s_nuc = []
for x in xlha:
    F2reps = []
    for i in range(1,251):
        apfel.SetReplica(i)
        apfel.ComputeStructureFunctionsAPFEL(Q0,Q)
        F2reps.append(apfel.F2total(x))
    F2s_nuc.append(F2reps)
        
F2s_nuc = np.array(F2s_nuc)    
errs_nuc = np.std(F2s_nuc, axis=1)
cvs_nuc = np.mean(F2s_nuc, axis=1)
    
cfac_deut = cvs_deuteron/cvs_proton
cfac_nuc = cvs_nuc/cvs_proton
errs_deut = errs_deuteron/cvs_proton
errs_nuc = errs_nuc/cvs_proton


# MMHT model - generate cfac
def generate_MMHT_cfac(xvals, N, c1, c2, c3e8, xp,
                       Nerr, c1err,  c2err, c3e8err):
    c1 = np.full_like(xvals, c1)
    c2 = np.full_like(xvals, c2)
    c3 = np.full_like(xvals, c3e8/1e8)
    xp = np.full_like(xvals, xp)
    c_4params = (1 + 0.01*N)*(1 + 0.01*c2*(np.log(xvals/xp))**2
                              + 0.01*c3*(np.log(xvals/xp))**20)
   
    # Calc error
    c2err = np.full_like(xvals, c2err)
    c3err = np.full_like(xvals, c3e8err/1e8)
    Nerr = np.full_like(xvals, Nerr)
    var = ((0.01*(1 + 0.01*c2*(np.log(xvals/xp))**2 +
                 0.01*c3*(np.log(xvals/xp))**20*Nerr))**2
          + (1 + 0.01*N)**2*((0.01*(np.log(xvals/xp))**2*c2err)**2
          + (0.01*(np.log(xvals/xp))**20*c3err)**2))
    err = np.sqrt(var)

    return c_4params, err

c_4params, err = generate_MMHT_cfac(xlha, 0.589, -0.116, -0.384, 0.0489, 0.03,
                                    0.738, 0.996, 0.182, 0.0056)


fig, ax = plt.subplots(figsize=(8,6))

plt.plot(xlha, cfac_deut, label="from fitted deuteron PDF", color="#756bb1", linewidth="5", linestyle=":")
plt.plot(xlha, cfac_nuc, label="from nNNPDF2.0", color="#2c7fb8", linestyle="--", linewidth="5")
plt.fill_between(xlha, cfac_deut-errs_deut, cfac_deut+errs_deut, alpha=0.5, color="#756bb1")
plt.fill_between(xlha,  cfac_nuc-errs_nuc, cfac_nuc+errs_nuc, alpha=0.5, color="#2c7fb8")
plt.plot(xlha, c_4params, "-", label="MMHT2014 NNLO 4 params", color="#a1d99b", linestyle="-", linewidth="5")
plt.fill_between(xlha, c_4params-err, c_4params+err, alpha=0.5, color="#a1d99b")
ax.set_xlim(0.01,1)
ax.set_xscale("log")
#ax.set_ylim([0.94,1.04])
xmin, xmax = ax.get_xlim()
ax.set_xlabel("x", fontsize="20")
ax.set_ylabel("Correction factor", fontsize="20")
plt.title(r"$Q^2$ = 10 GeV$^2$", fontsize="30")
ax.hlines(1, xmin, xmax, linestyles="-")
ax.legend(fontsize="15")



