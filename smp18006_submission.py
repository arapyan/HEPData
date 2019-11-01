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

sig_digits = 3
sig_digits2 = 2


submission.read_abstract("HEPData/inputs/smp18006/abstract.txt")
submission.add_link("Webpage with all figures and tables", "http://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-18-006/index.html")
submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:1905.07445")
submission.add_record_id(1735737, "inspire")

### Begin Table 2
table2 = Table("Table 2")
table2.description = "Expected yields from various background processes in $\mathrm{WV}$ and $\mathrm{ZV}$ final states. The combination of the statistical and systematic uncertainties are shown. The predicted yields are shown with their best-fit normalizations from the background-only fit. The aQGC signal yields are shown for two aQGC scenarios with $f_{T2}/ \Lambda^{4} = -0.5$ TeV$^{-4}$ and $f_{T2}/ \Lambda^{4} = -2.5$ TeV$^{-4}$ for the  $\mathrm{WV}$ and $\mathrm{ZV}$ channels, respectively. The charged Higgs boson signal yields are also shown for values of $s_{\mathrm{H}}=0.5$ and $m_{\mathrm{H}_{5}}=500$ GeV in the GM model. The statistical uncertainties are shown for the expected signal yields."
table2.location = "Data from Table 2"

table2.keywords["observables"] = ["Events"]

data2 = np.loadtxt("HEPData/inputs/smp18006/total_yields.txt", dtype='string', skiprows=2)

print(data2)

table2_data = Variable("Process", is_independent=True, is_binned=False, units="")
table2_data.values = [str(x) for x in data2[:,0]]

table2_yields1 = Variable("WV events", is_independent=False, is_binned=False, units="")
table2_yields1.values = [float(x) for x in data2[:,1]]
table2_yields1.add_qualifier("Expected events", "WV selection")
table2_yields1.add_qualifier("SQRT(S)", 13, "TeV")

table2_unc1 = Uncertainty("", is_symmetric=True)
table2_unc1.values = [float(x) for x in data2[:,2]]

table2_yields2 = Variable("ZV events", is_independent=False, is_binned=False, units="")
table2_yields2.values = [float(x) for x in data2[:,3]]
table2_yields2.add_qualifier("Expected events", "ZV selection")
table2_yields2.add_qualifier("SQRT(S)", 13, "TeV")

table2_unc2 = Uncertainty("", is_symmetric=True)
table2_unc2.values = [float(x) for x in data2[:,4]]

table2_yields1.add_uncertainty(table2_unc1)
table2_yields2.add_uncertainty(table2_unc2)


table2.add_variable(table2_data)
table2.add_variable(table2_yields1)
table2.add_variable(table2_yields2)

submission.add_table(table2)

for table2 in submission.tables:
    table2.keywords["cmenergies"] = [13000]

### End Table 2

### Begin Table 3
table3 = Table("Table 3")
table3.description = "Observed and expected lower and upper 95\% CL limits on the parameters of the quartic operators S0, S1, M0, M1, M6, M7, T0, T1, and T2 in $\mathrm{WV}$ and $\mathrm{ZV}$ channels. The last two columns show the observed and expected limits for the combination of the $\mathrm{WV}$ and $\mathrm{ZV}$ channels."
table3.location = "Data from Table 3"

table3.keywords["observables"] = ["Limits"]

data3 = np.loadtxt("HEPData/inputs/smp18006/limits.txt", dtype='string', skiprows=0)

print(data3)

table3_data = Variable("", is_independent=True, is_binned=False, units="")
table3_data.values = [str(x) for x in data3[:,0]]

table3_limits1 = Variable("Observed ($\mathrm{WV}$)", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits1.values = [str(x) for x in data3[:,1]]
#table3_limits1.add_qualifier("Expected events", "WV selection")
#table3_limits1.add_qualifier("SQRT(S)", 13, "TeV")

table3_limits2 = Variable("Expected ($\mathrm{WV}$)", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits2.values = [str(x) for x in data3[:,2]]

table3_limits3 = Variable("Observed ($\mathrm{ZV}$)", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits3.values = [str(x) for x in data3[:,3]]

table3_limits4 = Variable("Expected ($\mathrm{ZV}$)", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits4.values = [str(x) for x in data3[:,4]]

table3_limits5 = Variable("Observed", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits5.values = [str(x) for x in data3[:,5]]

table3_limits6 = Variable("Expected", is_independent=False, is_binned=False, units="TeV$^{-4}$")
table3_limits6.values = [str(x) for x in data3[:,6]]

table3.add_variable(table3_data)
table3.add_variable(table3_limits1)
table3.add_variable(table3_limits2)
table3.add_variable(table3_limits3)
table3.add_variable(table3_limits4)
table3.add_variable(table3_limits5)
table3.add_variable(table3_limits6)

submission.add_table(table3)

for table3 in submission.tables:
    table3.keywords["cmenergies"] = [13000]


### Begin Table limits
table_limits = Table("Table limits")
table_limits.description = "Expected and observed exclusion limits at the 95\% CL as a function of $m(\mathrm{H}^{\pm\pm})$ for $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm\pm}) \, \mathcal{B}(\mathrm{H}^{\pm\pm} \rightarrow \mathrm{W}^{\pm}\mathrm{W}^{\pm})$ in the $\mathrm{WV}$ final state."

table_limits.location = "Data from Figure 5-c"

table_limits.keywords["observables"] = ["Cross section in pb"]

data_limits = np.loadtxt("HEPData/inputs/smp18006/cross_section_limits.txt", dtype='float')

print(data_limits)

mass_limits = Variable("Mass", is_independent=True, is_binned=False, units="")
mass_limits.values = data_limits[0,:]

obs_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
obs_limits.values = data_limits[1,:]
obs_limits.add_qualifier("Cross section limits", "Observed limits")
obs_limits.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits.digits = sig_digits

exp_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
exp_limits.values = data_limits[2,:]
exp_limits.add_qualifier("Cross section limits", "Expected limits")
exp_limits.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits.digits = sig_digits

plus1sigma_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus1sigma_limits.values = data_limits[3,:]
plus1sigma_limits.add_qualifier("Cross section limits", "+1 sigma limits")
plus1sigma_limits.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits.digits = sig_digits

minus1sigma_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus1sigma_limits.values = data_limits[4,:]
minus1sigma_limits.add_qualifier("Cross section limits", "-1 sigma limits")
minus1sigma_limits.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits.digits = sig_digits

plus2sigma_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus2sigma_limits.values = data_limits[5,:]
plus2sigma_limits.add_qualifier("Cross section limits", "+2 sigma limits")
plus2sigma_limits.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits.digits = sig_digits

minus2sigma_limits = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus2sigma_limits.values = data_limits[6,:]
minus2sigma_limits.add_qualifier("Cross section limits", "-2 sigma limits")
minus2sigma_limits.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits.digits = sig_digits

table_limits.add_variable(mass_limits)
table_limits.add_variable(obs_limits)
table_limits.add_variable(exp_limits)
table_limits.add_variable(minus2sigma_limits)
table_limits.add_variable(minus1sigma_limits)
table_limits.add_variable(plus1sigma_limits)
table_limits.add_variable(plus2sigma_limits)

submission.add_table(table_limits)

for table_limits in submission.tables:
    table_limits.keywords["cmenergies"] = [13000]


### Begin Table limits 2
table_limits2 = Table("Table limits 2")
table_limits2.description = "Expected and observed exclusion limits at the 95\% CL as a function of $m(\mathrm{H}^{\pm})$ for $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm}) \, \mathcal{B}(\mathrm{H}^{\pm} \rightarrow \mathrm{W}^{\pm}\mathrm{Z})$ in the $\mathrm{ZV}$ final state."

table_limits2.location = "Data from Figure 5-b"

table_limits2.keywords["observables"] = ["Cross section in pb"]

data_limits2 = np.loadtxt("HEPData/inputs/smp18006/cross_section_limits2.txt", dtype='float')

print(data_limits2)

mass_limits2 = Variable("Mass", is_independent=True, is_binned=False, units="")
mass_limits2.values = data_limits2[0,:]

obs_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
obs_limits2.values = data_limits2[1,:]
obs_limits2.add_qualifier("Cross section limits", "Observed limits")
obs_limits2.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits2.digits = sig_digits

exp_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
exp_limits2.values = data_limits2[2,:]
exp_limits2.add_qualifier("Cross section limits", "Expected limits")
exp_limits2.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits2.digits = sig_digits

plus1sigma_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus1sigma_limits2.values = data_limits2[3,:]
plus1sigma_limits2.add_qualifier("Cross section limits", "+1 sigma limits")
plus1sigma_limits2.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits2.digits = sig_digits

minus1sigma_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus1sigma_limits2.values = data_limits2[4,:]
minus1sigma_limits2.add_qualifier("Cross section limits", "-1 sigma limits")
minus1sigma_limits2.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits2.digits = sig_digits

plus2sigma_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus2sigma_limits2.values = data_limits2[5,:]
plus2sigma_limits2.add_qualifier("Cross section limits", "+2 sigma limits")
plus2sigma_limits2.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits2.digits = sig_digits

minus2sigma_limits2 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus2sigma_limits2.values = data_limits2[6,:]
minus2sigma_limits2.add_qualifier("Cross section limits", "-2 sigma limits")
minus2sigma_limits2.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits2.digits = sig_digits

table_limits2.add_variable(mass_limits2)
table_limits2.add_variable(obs_limits2)
table_limits2.add_variable(exp_limits2)
table_limits2.add_variable(minus2sigma_limits2)
table_limits2.add_variable(minus1sigma_limits2)
table_limits2.add_variable(plus1sigma_limits2)
table_limits2.add_variable(plus2sigma_limits2)

submission.add_table(table_limits2)

for table_limits2 in submission.tables:
    table_limits2.keywords["cmenergies"] = [13000]


### Begin Table limits 3
table_limits3 = Table("Table limits 3")
table_limits3.description = "Expected and observed exclusion limits at the 95\% CL as a function of $m(\mathrm{H}^{\pm})$ for $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm}) \, \mathcal{B}(\mathrm{H}^{\pm} \rightarrow \mathrm{W}^{\pm}\mathrm{Z})$ in the $\mathrm{WV}$ final state."

table_limits3.location = "Data from Figure 5-a"

table_limits3.keywords["observables"] = ["Cross section in pb"]

data_limits3 = np.loadtxt("HEPData/inputs/smp18006/cross_section_limits3.txt", dtype='float')

print(data_limits3)

mass_limits3 = Variable("Mass", is_independent=True, is_binned=False, units="")
mass_limits3.values = data_limits3[0,:]

obs_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
obs_limits3.values = data_limits3[1,:]
obs_limits3.add_qualifier("Cross section limits", "Observed limits")
obs_limits3.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits3.digits = sig_digits

exp_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
exp_limits3.values = data_limits3[2,:]
exp_limits3.add_qualifier("Cross section limits", "Expected limits")
exp_limits3.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits3.diits = sig_digits

plus1sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus1sigma_limits3.values = data_limits3[3,:]
plus1sigma_limits3.add_qualifier("Cross section limits", "+1 sigma limits")
plus1sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits3.digits = sig_digits

minus1sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus1sigma_limits3.values = data_limits3[4,:]
minus1sigma_limits3.add_qualifier("Cross section limits", "-1 sigma limits")
minus1sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits3.digits = sig_digits

plus2sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus2sigma_limits3.values = data_limits3[5,:]
plus2sigma_limits3.add_qualifier("Cross section limits", "+2 sigma limits")
plus2sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits3.digits = sig_digits

minus2sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus2sigma_limits3.values = data_limits3[6,:]
minus2sigma_limits3.add_qualifier("Cross section limits", "-2 sigma limits")
minus2sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits3.digits = sig_digits


table_limits3.add_variable(mass_limits3)
table_limits3.add_variable(obs_limits3)
table_limits3.add_variable(exp_limits3)
table_limits3.add_variable(minus2sigma_limits3)
table_limits3.add_variable(minus1sigma_limits3)
table_limits3.add_variable(plus1sigma_limits3)
table_limits3.add_variable(plus2sigma_limits3)

submission.add_table(table_limits3)

for table_limits3 in submission.tables:
    table_limits3.keywords["cmenergies"] = [13000]


### Begin Table limits 4
table_limits4 = Table("Table limits 4")
table_limits4.description = "Expected and observed exclusion limits at the 95\% CL for $s_{\mathrm{H}}$ in the Georgi-Machacek model."

table_limits4.location = "Data from Figure 5-d"

table_limits4.keywords["observables"] = ["$s_{\mathrm{H}}$"]

data_limits4 = np.loadtxt("HEPData/inputs/smp18006/cross_section_limits_model.txt", dtype='float')

print(data_limits4)

mass_limits4 = Variable("Mass", is_independent=True, is_binned=False, units="")
mass_limits4.values = data_limits4[0,:]

obs_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
obs_limits4.values = data_limits4[1,:]
obs_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "Observed limits")
obs_limits4.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits4.digits = sig_digits2

exp_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
exp_limits4.values = data_limits4[2,:]
exp_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected limits")
exp_limits4.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits4.digits = sig_digits2

plus1sigma_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
plus1sigma_limits4.values = data_limits4[3,:]
plus1sigma_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "+1 sigma limits")
plus1sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits4.digits = sig_digits2

minus1sigma_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
minus1sigma_limits4.values = data_limits4[4,:]
minus1sigma_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "-1 sigma limits")
minus1sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits4.digits = sig_digits2

plus2sigma_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
plus2sigma_limits4.values = data_limits4[5,:]
plus2sigma_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "+2 sigma limits")
plus2sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits4.digits = sig_digits2

minus2sigma_limits4 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
minus2sigma_limits4.values = data_limits4[6,:]
minus2sigma_limits4.add_qualifier("$s_{\mathrm{H}}$ limits", "-2 sigma limits")
minus2sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits4.digits = sig_digits2

table_limits4.add_variable(mass_limits4)
table_limits4.add_variable(obs_limits4)
table_limits4.add_variable(exp_limits4)
table_limits4.add_variable(minus2sigma_limits4)
table_limits4.add_variable(minus1sigma_limits4)
table_limits4.add_variable(plus1sigma_limits4)
table_limits4.add_variable(plus2sigma_limits4)

submission.add_table(table_limits4)

for table_limits4 in submission.tables:
    table_limits4.keywords["cmenergies"] = [13000]

### Begin WV histograms
reader_wv = RootFileReader("HEPData/inputs/smp18006/wv_hist.root")

tableWV = Table("Figure 4-a")
tableWV.description = "$\mathrm{M_{\mathrm{WV}}}$ distribution in the $\mathrm{WV}$ channel. The predicted yields are shown with their best-fit normalizations from the background-only fit."
tableWV.location = "Data from Figure 4-a"
tableWV.keywords["observables"] = ["Events"]

histo_data_wv  = reader_wv.read_hist_1d("hdata")
histo_vvjjqcd_wv  = reader_wv.read_hist_1d("hVVjjQCD")
histo_vjet_wv  = reader_wv.read_hist_1d("hwjet")
histo_vvjj_wv  = reader_wv.read_hist_1d("hdiboson")
histo_top_wv  = reader_wv.read_hist_1d("htop")
histo_bkg_wv  = reader_wv.read_hist_1d("hbkg")

histo_data_wv.keys()

mmed_wv = Variable("$\mathrm{M_{\mathrm{WV}}}$", is_independent=True, is_binned=False, units="GeV")
mmed_wv.values = histo_vjet_wv["x"]

# y-axis: N events
data_wv = Variable("Number of data events", is_independent=False, is_binned=False, units="")
data_wv.values = histo_data_wv["y"]
data_wv.digits = sig_digits

#unc_data_wv = Uncertainty("Poisson errors", is_symmetric=True)
#unc_data_wv.values = histo_data_wv["dy"]
#data_wv.add_uncertainty(unc_data_wv)

vvjjqcd_wv = Variable("Number of QCD VV events", is_independent=False, is_binned=False, units="")
vvjjqcd_wv.values = histo_vvjjqcd_wv["y"]
vvjjqcd_wv.digits = sig_digits

vjet_wv = Variable("Number of V+jets events", is_independent=False, is_binned=False, units="")
vjet_wv.values = histo_vjet_wv["y"]
vjet_wv.digits = sig_digits

vvjj_wv = Variable("Number of SM EV VV events", is_independent=False, is_binned=False, units="")
vvjj_wv.values = histo_vvjj_wv["y"]
vvjj_wv.digits = sig_digits

top_wv = Variable("Number of top quark events", is_independent=False, is_binned=False, units="")
top_wv.values = histo_top_wv["y"]
top_wv.digits = sig_digits

background_wv = Variable("Number of background events", is_independent=False, is_binned=False, units="")
background_wv.values = histo_bkg_wv["y"]
background_wv.digits = sig_digits

unc_background_wv = Uncertainty("total uncertainty", is_symmetric=True)
unc_background_wv.values = histo_bkg_wv["dy"]
unc_background_wv.digits = 1
background_wv.add_uncertainty(unc_background_wv)

tableWV.add_variable(mmed_wv)
tableWV.add_variable(data_wv)
tableWV.add_variable(vvjjqcd_wv)
tableWV.add_variable(vjet_wv)
tableWV.add_variable(vvjj_wv)
tableWV.add_variable(top_wv)
tableWV.add_variable(background_wv)
submission.add_table(tableWV)
### End WV

### Begin ZV histograms
reader_zv = RootFileReader("HEPData/inputs/smp18006/zv_hist.root")

tableZV = Table("Figure 4-b")
tableZV.description = "$\mathrm{M_{\mathrm{ZV}}}$ distribution in the $\mathrm{ZV}$ channel. The predicted yields are shown with their best-fit normalizations from the background-only fit."
tableZV.location = "Data from Figure 4-b"
tableZV.keywords["observables"] = ["Events"]

histo_data_zv  = reader_zv.read_hist_1d("hdata")
histo_vvjjqcd_zv  = reader_zv.read_hist_1d("hVVjjQCD")
histo_vjet_zv  = reader_zv.read_hist_1d("hwjet")
histo_vvjj_zv  = reader_zv.read_hist_1d("hdiboson")
histo_top_zv  = reader_zv.read_hist_1d("htop")
histo_bkg_zv  = reader_zv.read_hist_1d("hbkg")

histo_data_zv.keys()

mmed_zv = Variable("$\mathrm{M_{\mathrm{ZV}}}$", is_independent=True, is_binned=False, units="GeV")
mmed_zv.values = histo_vjet_wv["x"]

# y-axis: N events
data_zv = Variable("Number of data events", is_independent=False, is_binned=False, units="")
data_zv.values = histo_data_zv["y"]
data_zv.digits = sig_digits

#unc_data_zv = Uncertainty("Poisson errors", is_symmetric=True)
#unc_data_zv.values = histo_data_zv["dy"]
#data_zv.add_uncertainty(unc_data_zv)

vvjjqcd_zv = Variable("Number of QCD VV events", is_independent=False, is_binned=False, units="")
vvjjqcd_zv.values = histo_vvjjqcd_zv["y"]
vvjjqcd_zv.digits = sig_digits

vjet_zv = Variable("Number of V+jets events", is_independent=False, is_binned=False, units="")
vjet_zv.values = histo_vjet_zv["y"]
vjet_zv.digits = sig_digits

vvjj_zv = Variable("Number of SM EV VV events", is_independent=False, is_binned=False, units="")
vvjj_zv.values = histo_vvjj_zv["y"]
vvjj_zv.digits = sig_digits

top_zv = Variable("Number of top quark events", is_independent=False, is_binned=False, units="")
top_zv.values = histo_top_zv["y"]
top_zv.digits = sig_digits

background_zv = Variable("Number of background events", is_independent=False, is_binned=False, units="")
background_zv.values = histo_bkg_zv["y"]
background_zv.digits = sig_digits

unc_background_zv = Uncertainty("total uncertainty", is_symmetric=True)
unc_background_zv.values = histo_bkg_zv["dy"]
unc_background_zv.digits = 1
background_zv.add_uncertainty(unc_background_zv)

tableZV.add_variable(mmed_zv)
tableZV.add_variable(data_zv)
tableZV.add_variable(vvjjqcd_zv)
tableZV.add_variable(vjet_zv)
tableZV.add_variable(vvjj_zv)
tableZV.add_variable(top_zv)
tableZV.add_variable(background_zv)
submission.add_table(tableZV)
### End ZV


### Begin covariance
# Create a reader for the input file
reader_covariance = RootFileReader("HEPData/inputs/smp18006/wv_mlfit.root")
# Read the histogram
data_covariance = reader_covariance.read_hist_2d("shapes_fit_b/ww13TeV2016/total_covar")
# Create variable objects
x_covariance = Variable("Bin X", is_independent=True, is_binned=False)
x_covariance.values = data_covariance["x"]
y_covariance = Variable("Bin Y", is_independent=True, is_binned=False)
y_covariance.values = data_covariance["y"]
z_covariance = Variable("covariance Matrix", is_independent=False, is_binned=False)
z_covariance.values = data_covariance["z"]

table_covariance = Table("Covariance Matrix")
table_covariance.description = "Covariance matrix for all the $\mathrm{M_{\mathrm{WV}}}$ bins used in the $\mathrm{WV}$ channel. The covariance matrix from the background-only fit is shown."
table_covariance.location = "Supplementary material"
for var in [x_covariance,y_covariance,z_covariance]:
    table_covariance.add_variable(var)
submission.add_table(table_covariance)
### End covariance


### Begin covariance ZV
# Create a reader for the input file
reader_covariance_zv = RootFileReader("HEPData/inputs/smp18006/zv_mlfit.root")
# Read the histogram
data_covariance_zv = reader_covariance_zv.read_hist_2d("shapes_fit_b/wzl13TeV2016/total_covar")
# Create variable objects
x_covariance_zv = Variable("Bin X", is_independent=True, is_binned=False)
x_covariance_zv.values = data_covariance["x"]
y_covariance_zv = Variable("Bin Y", is_independent=True, is_binned=False)
y_covariance_zv.values = data_covariance["y"]
z_covariance_zv = Variable("covariance Matrix", is_independent=False, is_binned=False)
z_covariance_zv.values = data_covariance_zv["z"]

table_covariance_zv = Table("Covariance Matrix 2")
table_covariance_zv.description = "Covariance matrix for all the $\mathrm{M_{\mathrm{ZV}}}$ bins used in the $\mathrm{ZV}$ channel. The covariance matrix from the background-only fit is shown."
table_covariance_zv.location = "Supplementary material"
for var in [x_covariance_zv,y_covariance_zv,z_covariance_zv]:
    table_covariance_zv.add_variable(var)
submission.add_table(table_covariance_zv)
### End covariance ZV


outdir = "smp18006_output"
submission.create_files(outdir)



