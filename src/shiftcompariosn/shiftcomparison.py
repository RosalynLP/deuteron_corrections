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
import lhapdf

xlha = np.logspace(-2, 0, num=100)[:-2]

# activate some options
apfel.SetPerturbativeOrder(2)
apfel.SetMassScheme("FONLL-C")
apfel.SetWMass(80.398)
apfel.SetProtonMass(0.938)
apfel.SetPoleMasses(1.51,4.92,172.5)
apfel.SetCKM(0.97427, 0.22536, 0.00355,0.2252, 0.97345, 0.041,0.00886, 0.04050, 0.99914)
apfel.EnableTargetMassCorrections(True)
apfel.EnableIntrinsicCharm(True)
apfel.SetTargetDIS("isoscalar")

# Setting proton only PDF set
apfel.SetPDFSet("NNPDF31_nnlo_as_0118_global_deut_ite")

# initializes integrals on the grids
apfel.InitializeAPFEL_DIS()

Q = 10

F2s_proton = []
for i in range(1,101):
    apfel.SetReplica(i)
    alphaQCD = lhapdf.mkAlphaS("NNPDF31_nnlo_as_0118_global_deut_ite",i)
    alphas = alphaQCD.alphasQ(Q)
    apfel.SetAlphaQCDRef(alphas,Q)
    apfel.ComputeStructureFunctionsAPFEL(Q,Q)
    F2reps = []
    for x in xlha:
        F2reps.append(apfel.F2total(x))
    F2s_proton.append(F2reps)
    
F2s_proton = np.array(F2s_proton)    
errs_proton = np.std(F2s_proton, axis=0)
cvs_proton = np.mean(F2s_proton, axis=0)

# Setting deuteron only PDF
apfel.SetPDFSet("NNPDF31_nnlo_as_0118_deuteron_only_ite")

F2s_deuteron = []
for i in range(1,101):
    apfel.SetReplica(i)
    alphaQCD = lhapdf.mkAlphaS("NNPDF31_nnlo_as_0118_deuteron_only_ite",i)
    alphas = alphaQCD.alphasQ(Q)
    apfel.SetAlphaQCDRef(alphas,Q)
    apfel.ComputeStructureFunctionsAPFEL(Q,Q)
    F2reps = []
    for x in xlha:
        F2reps.append(apfel.F2total(x))
    F2s_deuteron.append(F2reps)
        
F2s_deuteron = np.array(F2s_deuteron)    
errs_deuteron = np.std(F2s_deuteron, axis=0)
cvs_deuteron = np.mean(F2s_deuteron, axis=0)

# Nuclear fit deuteron PDF

apfel.SetPerturbativeOrder(1)
apfel.SetMassScheme("FONLL-B")
apfel.SetPDFSet("nNNPDF20_nlo_as_0118_p_A2_Z1")
apfel.EnableIntrinsicCharm(False)
#apfel.InitializeAPFEL_DIS()

F2s_nuc = []
for i in range(1,1001):
    apfel.SetReplica(i)
    alphaQCD = lhapdf.mkAlphaS("nNNPDF20_nlo_as_0118_p_A2_Z1",i)
    alphas = alphaQCD.alphasQ(Q)
    apfel.SetAlphaQCDRef(alphas,Q)
    apfel.ComputeStructureFunctionsAPFEL(Q,Q)
    F2reps = []
    for x in xlha:
        F2reps.append(apfel.F2total(x))
    F2s_nuc.append(F2reps)
        
F2s_nuc = np.array(F2s_nuc)    
errs_nuc = np.std(F2s_nuc, axis=0)
cvs_nuc = np.mean(F2s_nuc, axis=0)
    
cfac_deut = cvs_deuteron/cvs_proton
cfac_nuc = cvs_nuc/cvs_proton
errs_deut = np.sqrt((errs_deuteron/cvs_deuteron)**2. + (errs_proton/cvs_proton)**2.)*cfac_deut
errs_nuc  = np.sqrt((errs_nuc/cvs_nuc)**2. + (errs_proton/cvs_proton)**2.)*cfac_nuc

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

plt.plot(xlha, cfac_deut, label="deuteron-ite2 (NNLO)", color="cornflowerblue", linewidth="5", linestyle=":")
plt.plot(xlha, cfac_nuc, label="nNNPDF2.0 (NLO)", color="darkblue", linestyle="--", linewidth="5")
plt.fill_between(xlha, cfac_deut-errs_deut, cfac_deut+errs_deut, alpha=0.5, color="cornflowerblue")
plt.fill_between(xlha,  cfac_nuc-errs_nuc, cfac_nuc+errs_nuc, alpha=0.5, color="darkblue")
plt.plot(xlha, c_4params, "-", label="MMHT2014 (NNLO 4 params.)", color="orange", linestyle="-", linewidth="5")
plt.fill_between(xlha, c_4params-err, c_4params+err, alpha=0.5, color="orange")
ax.set_xlim(0.01,1)
ax.set_xscale("log")
ax.set_ylim([0.85,1.15])
xmin, xmax = ax.get_xlim()
ax.set_xlabel("x", fontsize="20")
#ax.set_ylabel("Correction factor", fontsize="20")
plt.title(r"Correction factor (Q=10 GeV)", fontsize="28")
ax.tick_params(labelsize="16")
ax.hlines(1, xmin, xmax, linestyles="-")
ax.legend(fontsize="15")
plt.savefig(f"deut_{Q}.png")



