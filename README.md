# deuteron_corrections
Code for determining deuteron uncertainties/corrections to NNPDF parton distribution functions. See [this accompanying paper](https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-08826-7.pdf).


## Contents of src
- PDFs

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code for plotting PDFs with deuteron uncertainties/corrections

- calc_chi2s

    Jupyter notebook for calculating partial chi2s per dataset with nuclear uncertainties/corrections

- c factor

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code to produce PDFs for c-factor computation

- chi2

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code for making chi2 reports

- covmat

    code for computing deuteron covariance matrix

- distances

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code for making PDF distance plots

- duplots

    code for computing the d/u quark ratio and plotting

- observables

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code for producing observables needed to compute deuteron covariance matrix. Also code for calculating and    plotting ratios of these observables

- rosalyn_pdfplots

    validphys runcards for use with [nnpdf](github.com/NNPDF/nnpdf) code for plotting PDFs
    
- shiftcompariosn

    script for comparing shift due to nuclear uncertainties with shift from MMHT2014


## Workflow

- Run validphys runcards in src/observables. Note that this requires a working version of [nnpdf](github.com/NNPDF/nnpdf).
- Run src/observables/observables.py to plot observable ratios
- Run src/covmat/covmat_new.py to create NNPDF4.0 style deuteron covariance matrices
- These can be included in a fit as a theory covariance matrix. See [nnpdf documentation](https://docs.nnpdf.science/).
- Plot PDFs using PDFs and rosalyn_pdfplots
- Compute c factors using cfactor
- Calculate chi2s using calc_chi2 and chi2
- Calculate d/u ratio using src/duplots
- Calculate PDF distances using src/distances
