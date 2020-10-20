#!/usr/bin/python

from __future__ import print_function
from hepdata_lib import Variable, Uncertainty
from hepdata_lib import Uncertainty
from hepdata_lib import RootFileReader

import hepdata_lib
from hepdata_lib import Submission

from hepdata_lib import Table
from hepdata_lib import Variable

import numpy as np
submission = Submission()


submission.read_abstract("HEPData/inputs/smp20006/abstract.txt")
submission.add_link("Webpage with all figures and tables", "http://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-20-006/index.html")
submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:2009.09429")
submission.add_record_id(1818160, "inspire")


### Begin Table 4
table4 = Table("Table 4")
table4.description = "Systematic uncertainties of the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{L}}$ and $\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$, and $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{X}}$ and $\mathrm{W}^\pm_{\mathrm{T}}\mathrm{W}^\pm_{\mathrm{T}}$ cross section measurements in units of percent."
table4.location = "Data from Table 4"

table4.keywords["observables"] = ["Uncertainty"]
table4.keywords["reactions"] = ["P P --> W W j j"]
table4.keywords["phrases"] = ["VBS", "Polarized", "Same-sign WW"]

data4 = np.loadtxt("HEPData/inputs/smp20006/systematics.txt", dtype='string', skiprows=2)

print(data4)

table4_data = Variable("Source of uncertainty", is_independent=True, is_binned=False, units="")
table4_data.values = [str(x) for x in data4[:,0]]

table4_yields0 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table4_yields0.values = [float(x) for x in data4[:,1]]
table4_yields0.add_qualifier("Source of uncertainty", "$\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{L}}$")
table4_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table4_yields0.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table4_yields1 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table4_yields1.values = [float(x) for x in data4[:,2]]
table4_yields1.add_qualifier("Source of uncertainty", "$\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$")
table4_yields1.add_qualifier("SQRT(S)", 13, "TeV")
table4_yields1.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table4_yields2 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table4_yields2.values = [float(x) for x in data4[:,3]]
table4_yields2.add_qualifier("Source of uncertainty", "$\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{X}}$")
table4_yields2.add_qualifier("SQRT(S)", 13, "TeV")
table4_yields2.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table4_yields3 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table4_yields3.values = [float(x) for x in data4[:,4]]
table4_yields3.add_qualifier("Source of uncertainty", "$\mathrm{W}^\pm_{\mathrm{T}}\mathrm{W}^\pm_{\mathrm{T}}$")
table4_yields3.add_qualifier("SQRT(S)", 13, "TeV")
table4_yields3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table4.add_variable(table4_data)
table4.add_variable(table4_yields0)
table4.add_variable(table4_yields1)
table4.add_variable(table4_yields2)
table4.add_variable(table4_yields3)

submission.add_table(table4)

for table4 in submission.tables:
    table4.keywords["cmenergies"] = [13000]
### End Table 4


### Begin Table 5
table5 = Table("Table 5")
table5.description = "Expected yields from various SM processes and observed data events in WW SR. The combination of the statistical and systematic uncertainties is shown. The expected yields are shown with their best fit normalizations from the simultaneous fit for the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{L}}$ and $\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$ cross sections. The $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{T}}$ and $\mathrm{W}^\pm_{\mathrm{T}}\mathrm{W}^\pm_{\mathrm{T}}$ yields are obtained from the $\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$ yield assuming the SM prediction for the ratio of the yields. The tVx background yield includes the contributions from  tt$\mathrm{V}$ and tZq processes."

table5.location = "Data from Table 5"

table5.keywords["observables"] = ["Events"]
table5.keywords["reactions"] = ["P P --> W W j j"]
table5.keywords["phrases"] = ["VBS", "Polarized", "Same-sign WW"]

data5 = np.loadtxt("HEPData/inputs/smp20006/total_yields_llfit.txt", dtype='string')

print(data5)

table5_data = Variable("Process", is_independent=True, is_binned=False, units="")
table5_data.values = [str(x) for x in data5[:,0]]

table5_yields0 = Variable("Yields in WW signal region", is_independent=False, is_binned=False, units="")
table5_yields0.values = [float(x) for x in data5[:,1]]
table5_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table5_unc0 = Uncertainty("total uncertainty", is_symmetric=True)
table5_unc0.values = [float(x) for x in data5[:,2]]
table5_yields0.add_uncertainty(table5_unc0)
table5_yields0.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table5.add_variable(table5_data)
table5.add_variable(table5_yields0)

submission.add_table(table5)

for table5 in submission.tables:
    table5.keywords["cmenergies"] = [13000]
###End Table 5

### Begin Table 6
table6 = Table("Table 6")
table6.description = "Measured fiducial cross sections for the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{L}}$ and $\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$  processes, and for the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{X}}$ and $\mathrm{W}^\pm_{\mathrm{T}}\mathrm{W}^\pm_{\mathrm{T}}$ processes for the helicity eigenstates defined in the WW center-of-mass frame. The combination of the statistical and systematic uncertainties is shown. The theoretical uncertainties include statistical, PDF, and LO scale uncertainties.  $\mathcal{B}$ is the branching fraction for $\mathrm{W}\mathrm{W} \\rightarrow \ell \\nu \ell' \\nu$. The fiducial region is defined by requiring two same-sign leptons with $p_{T}>20$, $|\eta|<2.5$, and $m_{ll}>20$, and two jets with $m_{jj}>500$ and $|\Delta \eta_{jj}|>2.5$. The jets at generator level are clustered from stable particles, excluding neutrinos, using the anti-kt clustering algorithm with R = 0.4, and are required to have $p_{T}>50$ and $|\eta|<4.7$. The jets within $\Delta R<0.4$ of the selected charged leptons are not included."
table6.location = "Data from Table 6"

table6.keywords["observables"] = ["Fiducial cross section"]
table6.keywords["reactions"] = ["P P --> W W j j"]
table6.keywords["phrases"] = ["VBS", "Polarized", "Same-sign WW"]

data6 = np.loadtxt("HEPData/inputs/smp20006/cross_sections_ww.txt", dtype='string')

print(data6)

table6_data = Variable("Process", is_independent=True, is_binned=False, units="")
table6_data.values = [str(x) for x in data6[:,0]]

table6_yields0 = Variable("$\sigma\mathcal{B}$", is_independent=False, is_binned=False, units="fb")
table6_yields0.values = [float(x) for x in data6[:,1]]
table6_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table6_unc0 = Uncertainty("total uncertainty", is_symmetric=False)
table6_unc0.values = [eval(x) for x in data6[:,2]]
table6_yields0.add_uncertainty(table6_unc0)

table6_yields1 = Variable("Theoretical prediction", is_independent=False, is_binned=False, units="fb")
table6_yields1.values = [float(x) for x in data6[:,3]]
table6_yields1.add_qualifier("SQRT(S)", 13, "TeV")
table6_unc1 = Uncertainty("total uncertainty", is_symmetric=True)
table6_unc1.values = [float(x) for x in data6[:,4]]
table6_yields1.add_uncertainty(table6_unc1)

table6.add_variable(table6_data)
table6.add_variable(table6_yields0)
table6.add_variable(table6_yields1)

submission.add_table(table6)

for table6 in submission.tables:
    table6.keywords["cmenergies"] = [13000]

### End Table 6

### Begin Table 7
table7 = Table("Table 7")
table7.description = "Measured fiducial cross sections for the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{L}}$ and $\mathrm{W}^\pm_{\mathrm{X}}\mathrm{W}^\pm_{\mathrm{T}}$  processes, and for the $\mathrm{W}^\pm_{\mathrm{L}}\mathrm{W}^\pm_{\mathrm{X}}$ and $\mathrm{W}^\pm_{\mathrm{T}}\mathrm{W}^\pm_{\mathrm{T}}$ processes for the helicity eigenstates defined in the parton-parton center-of-mass frame. The combination of the statistical and systematic uncertainties is shown. The theoretical uncertainties include statistical, PDF, and LO scale uncertainties.  $\mathcal{B}$ is the branching fraction for $\mathrm{W}\mathrm{W} \\rightarrow \ell \\nu \ell' \\nu$. The fiducial region is defined by requiring two same-sign leptons with $p_{T}>20$, $|\eta|<2.5$, and $m_{ll}>20$, and two jets with $m_{jj}>500$ and $|\Delta \eta_{jj}|>2.5$. The jets at generator level are clustered from stable particles, excluding neutrinos, using the anti-kt clustering algorithm with R = 0.4, and are required to have $p_{T}>50$ and $|\eta|<4.7$. The jets within $\Delta R<0.4$ of the selected charged leptons are not included."
table7.location = "Data from Table 7"

table7.keywords["observables"] = ["Fiducial cross section"]
table7.keywords["reactions"] = ["P P --> W W j j"]
table7.keywords["phrases"] = ["VBS", "Polarized", "Same-sign WW"]

data7 = np.loadtxt("HEPData/inputs/smp20006/cross_sections_pp.txt", dtype='string')

print(data7)

table7_data = Variable("Process", is_independent=True, is_binned=False, units="")
table7_data.values = [str(x) for x in data7[:,0]]

table7_yields0 = Variable("$\sigma\mathcal{B}$", is_independent=False, is_binned=False, units="fb")
table7_yields0.values = [float(x) for x in data7[:,1]]
table7_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table7_unc0 = Uncertainty("total uncertainty", is_symmetric=False)
table7_unc0.values = [eval(x) for x in data7[:,2]]
table7_yields0.add_uncertainty(table7_unc0)

table7_yields1 = Variable("Theoretical prediction", is_independent=False, is_binned=False, units="fb")
table7_yields1.values = [float(x) for x in data7[:,3]]
table7_yields1.add_qualifier("SQRT(S)", 13, "TeV")
table7_unc1 = Uncertainty("total uncertainty", is_symmetric=True)
table7_unc1.values = [float(x) for x in data7[:,4]]
table7_yields1.add_uncertainty(table7_unc1)

table7.add_variable(table7_data)
table7.add_variable(table7_yields0)
table7.add_variable(table7_yields1)

submission.add_table(table7)

for table7 in submission.tables:
    table7.keywords["cmenergies"] = [13000]

### End Table 7


#
outdir = "smp20006_output"
submission.create_files(outdir)



