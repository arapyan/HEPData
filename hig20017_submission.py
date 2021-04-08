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

submission.read_abstract("HEPData/inputs/hig20017/abstract.txt")
submission.add_link("Webpage with all figures and tables", "http://cms-results.web.cern.ch/cms-results/public-results/publications/HIG-20-017/index.html")
#submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:xxxx.xxxxx")
#submission.add_record_id(1818160, "inspire")


### Begin Table 2
table2 = Table("Table 2")
table2.description = "Summary of the impact of the systematic uncertainties on the extracted signal strength; for the case of a background-only simulated data set, i.e., assuming no contributions from the $\mathrm{H}^{\pm}$ and $\mathrm{H}^{\pm\pm}$ processes, and including a charged Higgs boson signal for values of $s_{\mathrm{H}}=1.0$ and $m_{\mathrm{H}_{5}}=500$ GeV in the GM model."
table2.location = "Data from Table 2"

table2.keywords["observables"] = ["Uncertainty"]
table2.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
table2.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data2 = np.loadtxt("HEPData/inputs/hig20017/systematics.txt", dtype='string', skiprows=2)

print(data2)

table2_data = Variable("Source of uncertainty", is_independent=True, is_binned=False, units="")
table2_data.values = [str(x) for x in data2[:,0]]

table2_yields0 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table2_yields0.values = [float(x) for x in data2[:,1]]
table2_yields0.add_qualifier("Source of uncertainty", "$\Delta \mu$ for background-only")
table2_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table2_yields0.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table2_yields1 = Variable("Uncertainty", is_independent=False, is_binned=False, units="")
table2_yields1.values = [float(x) for x in data2[:,2]]
table2_yields1.add_qualifier("Source of uncertainty", "$\Delta \mu$ for $s_{\mathrm{H}}=1.0$ and $m_{\mathrm{H}_{5}}=500~\mathrm{GeV}$")
table2_yields1.add_qualifier("SQRT(S)", 13, "TeV")
table2_yields1.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table2.add_variable(table2_data)
table2.add_variable(table2_yields0)
table2.add_variable(table2_yields1)

submission.add_table(table2)

for table2 in submission.tables:
    table2.keywords["cmenergies"] = [13000]

### End Table 2


### Begin Table 3
table3 = Table("Table 3")
table3.description = "Expected signal and background yields from various SM processes and observed data events in all regions used in the analysis. The expected background yields are shown with their normalizations from the simultaneous fit for the background-only hypothesis, i.e., assuming no contributions from the  $\mathrm{H}^{\pm}$  and $\mathrm{H}^{\pm\pm}$ processes. The expected signal yields are shown for $s_{\mathrm{H}}=1.0$ in the GM model. The combination of the statistical and systematic uncertainties is shown."
table3.location = "Data from Table 3"

table3.keywords["observables"] = ["Events"]
table3.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
table3.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data3 = np.loadtxt("HEPData/inputs/hig20017/yields.txt", dtype='string', skiprows=2)

print(data3)

table3_data = Variable("Process", is_independent=True, is_binned=False, units="")
table3_data.values = [str(x) for x in data3[:,0]]

table3_yields0 = Variable("$\mathrm{W}^{\pm}\mathrm{W}^{\pm}$ SR", is_independent=False, is_binned=False, units="")
table3_unc0 = Uncertainty("total uncertainty", is_symmetric=True)
table3_yields0.values = [float(x) for x in data3[:,1]]
table3_unc0.values = [float(x) for x in data3[:,2]]
table3_yields0.add_uncertainty(table3_unc0)
table3_yields0.add_qualifier("SQRT(S)", 13, "TeV")
table3_yields0.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table3_yields1 = Variable("WZ SR", is_independent=False, is_binned=False, units="")
table3_unc1 = Uncertainty("total uncertainty", is_symmetric=True)
table3_yields1.values = [float(x) for x in data3[:,3]]
table3_unc1.values = [float(x) for x in data3[:,4]]
table3_yields1.add_uncertainty(table3_unc1)
table3_yields1.add_qualifier("SQRT(S)", 13, "TeV")
table3_yields1.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table3_yields2 = Variable("Nonprompt CR", is_independent=False, is_binned=False, units="")
table3_unc2 = Uncertainty("total uncertainty", is_symmetric=True)
table3_yields2.values = [float(x) for x in data3[:,5]]
table3_unc2.values = [float(x) for x in data3[:,6]]
table3_yields2.add_uncertainty(table3_unc2)
table3_yields2.add_qualifier("SQRT(S)", 13, "TeV")
table3_yields2.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table3_yields3 = Variable("tZq CR", is_independent=False, is_binned=False, units="")
table3_unc3 = Uncertainty("total uncertainty", is_symmetric=True)
table3_yields3.values = [float(x) for x in data3[:,7]]
table3_unc3.values = [float(x) for x in data3[:,8]]
table3_yields3.add_uncertainty(table3_unc3)
table3_yields3.add_qualifier("SQRT(S)", 13, "TeV")
table3_yields3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table3_yields4 = Variable("ZZ CR", is_independent=False, is_binned=False, units="")
table3_unc4 = Uncertainty("total uncertainty", is_symmetric=True)
table3_yields4.values = [float(x) for x in data3[:,9]]
table3_unc4.values = [float(x) for x in data3[:,10]]
table3_yields4.add_uncertainty(table3_unc4)
table3_yields4.add_qualifier("SQRT(S)", 13, "TeV")
table3_yields4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")

table3.add_variable(table3_data)
table3.add_variable(table3_yields0)
table3.add_variable(table3_yields1)
table3.add_variable(table3_yields2)
table3.add_variable(table3_yields3)
table3.add_variable(table3_yields4)

submission.add_table(table3)

for table3 in submission.tables:
    table3.keywords["cmenergies"] = [13000]

### End Table 2

### Begin Fig4
reader_Fig4 = RootFileReader("HEPData/inputs/hig20017/histoDatacard_hpp_hp_fiducial6_mH500_2019.root")
reader_fit = RootFileReader("HEPData/inputs/hig20017/fitDiagnosticsssww_2019_fiducial6_mH500_obs.root")


tableFig4 = Table("Figure 4")
tableFig4.description = "Distributions for signal, backgrounds, and data for the bins used in the simultaneous fit. The bins 1--32 (4$\\times$8) show the events in the WW SR ($m_{\mathrm{jj}} \\times m_{\mathrm{T}}$), the bins 33--46 (2$\\times$7) show the events in the WZ SR ($m_{\mathrm{jj}} \\times m_{\mathrm{T}}$), the 4 bins 47--50 show the events in the nonprompt lepton CR ($m_{\mathrm{jj}}$), the 4 bins 51--54 show the events in the tZq CR ($m_{\mathrm{jj}}$), and the 4 bins 55--58 show the events in the ZZ CR ($m_{\mathrm{jj}}$). The predicted yields are shown with their best fit normalizations from the simultaneous fit  for the background-only hypothesis, i.e., assuming no contributions from the  $\mathrm{H}^{\pm}$ and $\mathrm{H}^{\pm\pm}$ processes. Vertical bars on data points represent the statistical uncertainty in the data. The histograms for tVx backgrounds include the contributions from ttV and tZq processes. The histograms for other backgrounds  include the contributions from double parton scattering, VVV, and from oppositely charged dilepton final states from tt, tW, $\mathrm{W}^{+}\mathrm{W}^{-}$, and Drell--Yan processes. The overflow is included in the last bin in each corresponding region. The lower panels show the ratio of the number of events observed in data to that of the total SM prediction. The hatched gray bands represent the uncertainties in the predicted yields. The solid lines show the signal predictions for values of $s_{\mathrm{H}}=1.0$ and $m_{\mathrm{H}_{5}}=500$ GeV in the GM model."
tableFig4.location = "Data from Figure 4"
tableFig4.keywords["observables"] = ["Events"]
tableFig4.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
tableFig4.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

histo0_Fig4    = reader_Fig4.read_hist_1d("histo0")  # data
histo5_Fig4    = reader_Fig4.read_hist_1d("histo5")  # WpWp
histo8_Fig4    = reader_Fig4.read_hist_1d("histo8")  # WZ
histo9_Fig4    = reader_Fig4.read_hist_1d("histo9")  # ZZ
histo10_Fig4   = reader_Fig4.read_hist_1d("histo10") # Nonprompt
histo12_Fig4   = reader_Fig4.read_hist_1d("histo12") # TVX
histo18_Fig4   = reader_Fig4.read_hist_1d("histo18") # Other
histo19_Fig4   = reader_Fig4.read_hist_1d("histo19") # H++
histo21_Fig4   = reader_Fig4.read_hist_1d("histo21") # H+
histo_Bck_Fig4 = reader_fit.read_hist_1d("shapes_fit_b/total_background")

histo0_Fig4.keys()

mmed_Fig4 = Variable("Bin", is_independent=True, is_binned=True, units="")
mmed_Fig4.values = histo0_Fig4["x_edges"]

totalbackground_Fig4 = Variable("Number of background events", is_independent=False, is_binned=False, units="")
totalbackground_Fig4.values = histo_Bck_Fig4["y"]

unc_totalbackground_Fig4 = Uncertainty("total uncertainty", is_symmetric=True)
unc_totalbackground_Fig4.values = histo_Bck_Fig4["dy"]

totalbackground_Fig4.add_uncertainty(unc_totalbackground_Fig4)

WW_Fig4 = Variable("Number of $\mathrm{W}^{\pm}\mathrm{W}^{\pm}$ events", is_independent=False, is_binned=False, units="")
WW_Fig4.values = histo5_Fig4["y"]

WZ_Fig4 = Variable("Number of WZ events", is_independent=False, is_binned=False, units="")
WZ_Fig4.values = histo8_Fig4["y"]

ZZ_Fig4 = Variable("Number of ZZ events", is_independent=False, is_binned=False, units="")
ZZ_Fig4.values = histo9_Fig4["y"]

Fake_Fig4 = Variable("Number of nonprompt lepton events", is_independent=False, is_binned=False, units="")
Fake_Fig4.values = histo10_Fig4["y"]

TVX_Fig4 = Variable("Number of tVx events", is_independent=False, is_binned=False, units="")
TVX_Fig4.values = histo12_Fig4["y"]

Other_Fig4 = Variable("Number of Other events", is_independent=False, is_binned=False, units="")
Other_Fig4.values = histo18_Fig4["y"]

Data_Fig4 = Variable("Number of data events", is_independent=False, is_binned=False, units="")
Data_Fig4.values = histo0_Fig4["y"]

hpp_Fig4 = Variable("Number of $\mathrm{H}^{\pm\pm}(500)\\rightarrow\mathrm{W}^{\pm}\mathrm{W}^{\pm}$ events", is_independent=False, is_binned=False, units="")
hpp_Fig4.values = histo19_Fig4["y"]

unc_hpp_Fig4 = Uncertainty("total uncertainty", is_symmetric=True)
unc_hpp_Fig4.values = histo19_Fig4["dy"]
hpp_Fig4.add_uncertainty(unc_hpp_Fig4)

hp_Fig4 = Variable("Number of $\mathrm{H}^{\pm}(500)\\rightarrow\mathrm{W}^{\pm}\mathrm{Z}$ events", is_independent=False, is_binned=False, units="")
hp_Fig4.values = histo21_Fig4["y"]

unc_hp_Fig4 = Uncertainty("total uncertainty", is_symmetric=True)
unc_hp_Fig4.values = histo21_Fig4["dy"]
hp_Fig4.add_uncertainty(unc_hp_Fig4)

#unc_data_Fig3a = Uncertainty("Poisson errors", is_symmetric=True)
#unc_data_Fig3a.values = histo0_Fig3a["dy"]
#Data_Fig3a.add_uncertainty(unc_data_Fig3a)

tableFig4.add_variable(mmed_Fig4)
tableFig4.add_variable(totalbackground_Fig4)
tableFig4.add_variable(WW_Fig4)
tableFig4.add_variable(WZ_Fig4)
tableFig4.add_variable(ZZ_Fig4)
tableFig4.add_variable(Fake_Fig4)
tableFig4.add_variable(TVX_Fig4)
tableFig4.add_variable(Other_Fig4)
tableFig4.add_variable(Data_Fig4)
tableFig4.add_variable(hpp_Fig4)
tableFig4.add_variable(hp_Fig4)
submission.add_table(tableFig4)
### End Fig4

### Begin Fig5
fig5 = Table("Figure 5")
fig5.description = "The product of acceptance and selection efficiency within the fiducial region for the VBF $\mathrm{H}^{\pm\pm}\\rightarrow\mathrm{W}^{\pm}\mathrm{W}^{\pm}\\rightarrow 2\ell 2v$ and $\mathrm{H}^{\pm}\\rightarrow\mathrm{W}^{\pm}\mathrm{Z}\\rightarrow 3\ell v$ processes, as a function of $m_{\mathrm{H}_{5}}$. The combination of the statistical and systematic uncertainties is shown. The theoretical uncertainties in the acceptance are also included."
fig5.location = "Data from Figure 5"

fig5.keywords["observables"] = ["Acceptance times efficiency"]
fig5.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
fig5.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data5 = np.loadtxt("HEPData/inputs/hig20017/hpp_eff.txt", dtype='string', skiprows=0)

print(data5)

fig5_data = Variable("$m_{\mathrm{H}_5}$", is_independent=True, is_binned=False, units="GeV")
fig5_data.values = [float(x) for x in data5[:,0]]

fig5_yields0 = Variable("Acceptance times efficiency for $\mathrm{H}^{\pm\pm}\\rightarrow\mathrm{W}^{\pm}\mathrm{W}^{\pm}\\rightarrow 2\ell 2v$", is_independent=False, is_binned=False, units="")
fig5_unc0 = Uncertainty("total uncertainty", is_symmetric=True)
fig5_yields0.values = [float(x) for x in data5[:,1]]
fig5_unc0.values = [eval(x) for x in data5[:,3]]
fig5_yields0.add_uncertainty(fig5_unc0)

fig5_yields1 = Variable("Acceptance times efficiency for $\mathrm{H}^{\pm}\\rightarrow\mathrm{W}^{\pm}\mathrm{Z}\\rightarrow 3\ell v$", is_independent=False, is_binned=False, units="")
fig5_unc1 = Uncertainty("total uncertainty", is_symmetric=True)
fig5_yields1.values = [float(x) for x in data5[:,2]]
fig5_unc1.values = [eval(x) for x in data5[:,4]]
fig5_yields1.add_uncertainty(fig5_unc1)

fig5.add_variable(fig5_data)
fig5.add_variable(fig5_yields0)
fig5.add_variable(fig5_yields1)

submission.add_table(fig5)

for fig5 in submission.tables:
    fig5.keywords["cmenergies"] = [13000]

### End Table 2

### Begin Table limits 3
table_limits3 = Table("Table limits 6-a")
table_limits3.description = "Expected and observed exclusion limits at 95\% confidence level for $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm\pm}) \mathcal{B}(\mathrm{H}^{\pm\pm} \\rightarrow \mathrm{W}^{\pm}\mathrm{W}^{\pm})$ as functions of $m_{\mathrm{H}^{\pm\pm}}$. The contribution of the $\mathrm{H}^{\pm}$ boson signal is set to zero for the derivation of the exclusion limits on the $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm\pm}) \mathcal{B}(\mathrm{H}^{\pm\pm} \\rightarrow \mathrm{W}^{\pm}\mathrm{W}^{\pm})$. Values above the curves are excluded."

table_limits3.location = "Data from Figure 6-a"

table_limits3.keywords["observables"] = ["Limits"]
table_limits3.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
table_limits3.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data_limits3 = np.loadtxt("HEPData/inputs/hig20017/hpp_limits.txt", dtype='float')

print(data_limits3)

mass_limits3 = Variable("Mass", is_independent=True, is_binned=False, units="GeV")
mass_limits3.values = data_limits3[:,0]

obs_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
obs_limits3.values = data_limits3[:,1]
obs_limits3.add_qualifier("Cross section limits", "Observed limits")
obs_limits3.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
obs_limits3.digits = sig_digits

exp_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
exp_limits3.values = data_limits3[:,3]
exp_limits3.add_qualifier("Cross section limits", "Expected limits")
exp_limits3.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
exp_limits3.digits = sig_digits

plus1sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus1sigma_limits3.values = data_limits3[:,6]
plus1sigma_limits3.add_qualifier("Cross section limits", "Expected +1 sigma limits")
plus1sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus1sigma_limits3.digits = sig_digits

minus1sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus1sigma_limits3.values = data_limits3[:,5]
minus1sigma_limits3.add_qualifier("Cross section limits", "Expected -1 sigma limits")
minus1sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
minus1sigma_limits3.digits = sig_digits

plus2sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus2sigma_limits3.values = data_limits3[:,7]
plus2sigma_limits3.add_qualifier("Cross section limits", "Expected +2 sigma limits")
plus2sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus2sigma_limits3.digits = sig_digits

minus2sigma_limits3 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus2sigma_limits3.values = data_limits3[:,4]
minus2sigma_limits3.add_qualifier("Cross section limits", "Expected -2 sigma limits")
minus2sigma_limits3.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits3.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
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
table_limits4 = Table("Table limits 6-b")
table_limits4.description = "Expected and observed exclusion limits at 95\% confidence level for $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm}) \mathcal{B}(\mathrm{H}^{\pm} \\rightarrow \mathrm{W}^{\pm}\mathrm{Z})$ as functions of $m_{\mathrm{H}^{\pm}}$. The contribution of the $\mathrm{H}^{\pm\pm}$ boson signal is set to zero for the derivation of the exclusion limits on the $\sigma_\mathrm{VBF}(\mathrm{H}^{\pm}) \mathcal{B}(\mathrm{H}^{\pm} \\rightarrow \mathrm{W}^{\pm}\mathrm{Z})$. Values above the curves are excluded."

table_limits4.location = "Data from Figure 6-b"

table_limits4.keywords["observables"] = ["Limits"]
table_limits4.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
table_limits4.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data_limits4 = np.loadtxt("HEPData/inputs/hig20017/hp_limits.txt", dtype='float')

print(data_limits4)

mass_limits4 = Variable("Mass", is_independent=True, is_binned=False, units="GeV")
mass_limits4.values = data_limits4[:,0]

obs_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
obs_limits4.values = data_limits4[:,1]
obs_limits4.add_qualifier("Cross section limits", "Observed limits")
obs_limits4.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
obs_limits4.digits = sig_digits

exp_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
exp_limits4.values = data_limits4[:,3]
exp_limits4.add_qualifier("Cross section limits", "Expected limits")
exp_limits4.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
exp_limits4.digits = sig_digits

plus1sigma_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus1sigma_limits4.values = data_limits4[:,6]
plus1sigma_limits4.add_qualifier("Cross section limits", "Expected +1 sigma limits")
plus1sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus1sigma_limits4.digits = sig_digits

minus1sigma_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus1sigma_limits4.values = data_limits4[:,5]
minus1sigma_limits4.add_qualifier("Cross section limits", "Expected -1 sigma limits")
minus1sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
minus1sigma_limits4.digits = sig_digits

plus2sigma_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
plus2sigma_limits4.values = data_limits4[:,7]
plus2sigma_limits4.add_qualifier("Cross section limits", "Expected +2 sigma limits")
plus2sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus2sigma_limits4.digits = sig_digits

minus2sigma_limits4 = Variable("Cross section", is_independent=False, is_binned=False, units="pb")
minus2sigma_limits4.values = data_limits4[:,4]
minus2sigma_limits4.add_qualifier("Cross section limits", "Expected -2 sigma limits")
minus2sigma_limits4.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits4.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
minus2sigma_limits4.digits = sig_digits

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

### Begin Table limits 5
table_limits5 = Table("Table limits 6-c")
table_limits5.description = "Expected and observed exclusion limits at 95\% confidence level for $s_{\mathrm{H}}$ as functions of $m_{\mathrm{H}_{5}}$ in the Georgi-Machacek model. Values above the curves are excluded. The exclusion limits for $s_{\mathrm{H}}$ are shown up to $m_{\mathrm{H}_{5}}=2000$ GeV, given the low sensitivity in the Georgi-Machacek model for values above that mass."

table_limits5.location = "Data from Figure 6-c"

table_limits5.keywords["observables"] = ["Limits"]
table_limits5.keywords["reactions"] = ["P P --> W W j j", "P P --> W Z j j"]
table_limits5.keywords["phrases"] = ["Same-sign WW", "WZ", "Georgi-Machacek", "Charged Higgs", "VBF"]

data_limits5 = np.loadtxt("HEPData/inputs/hig20017/gm_limits.txt", dtype='float')

print(data_limits5)

mass_limits5 = Variable("Mass", is_independent=True, is_binned=False, units="GeV")
mass_limits5.values = data_limits5[:,0]

obs_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
obs_limits5.values = data_limits5[:,1]
obs_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Observed limits")
obs_limits5.add_qualifier("SQRT(S)", 13, "TeV")
obs_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
obs_limits5.digits = sig_digits

exp_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
exp_limits5.values = data_limits5[:,3]
exp_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected limits")
exp_limits5.add_qualifier("SQRT(S)", 13, "TeV")
exp_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
exp_limits5.digits = sig_digits

plus1sigma_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
plus1sigma_limits5.values = data_limits5[:,6]
plus1sigma_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected +1 sigma limits")
plus1sigma_limits5.add_qualifier("SQRT(S)", 13, "TeV")
plus1sigma_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus1sigma_limits5.digits = sig_digits

minus1sigma_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
minus1sigma_limits5.values = data_limits5[:,5]
minus1sigma_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected -1 sigma limits")
minus1sigma_limits5.add_qualifier("SQRT(S)", 13, "TeV")
minus1sigma_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
minus1sigma_limits5.digits = sig_digits

plus2sigma_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
plus2sigma_limits5.values = data_limits5[:,7]
plus2sigma_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected +2 sigma limits")
plus2sigma_limits5.add_qualifier("SQRT(S)", 13, "TeV")
plus2sigma_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
plus2sigma_limits5.digits = sig_digits

minus2sigma_limits5 = Variable("$s_{\mathrm{H}}$", is_independent=False, is_binned=False, units="")
minus2sigma_limits5.values = data_limits5[:,4]
minus2sigma_limits5.add_qualifier("$s_{\mathrm{H}}$ limits", "Expected -2 sigma limits")
minus2sigma_limits5.add_qualifier("SQRT(S)", 13, "TeV")
minus2sigma_limits5.add_qualifier("L$_{\mathrm{int}}$", 137, "fb$^{-1}$")
minus2sigma_limits5.digits = sig_digits

table_limits5.add_variable(mass_limits5)
table_limits5.add_variable(obs_limits5)
table_limits5.add_variable(exp_limits5)
table_limits5.add_variable(minus2sigma_limits5)
table_limits5.add_variable(minus1sigma_limits5)
table_limits5.add_variable(plus1sigma_limits5)
table_limits5.add_variable(plus2sigma_limits5)

submission.add_table(table_limits5)

for table_limits5 in submission.tables:
    table_limits5.keywords["cmenergies"] = [13000]


### Begin covariance

# Read the histogram
data_covariance = reader_fit.read_hist_2d("shapes_fit_b/overall_total_covar")
# Create variable objects
x_covariance = Variable("Bin X", is_independent=True, is_binned=True)
x_covariance.values = data_covariance["x_edges"]
y_covariance = Variable("Bin Y", is_independent=True, is_binned=True)
y_covariance.values = data_covariance["y_edges"]
z_covariance = Variable("covariance matrix", is_independent=False, is_binned=False)
z_covariance.values = data_covariance["z"]

table_covariance = Table("Covariance Matrix")
table_covariance.description = "Covariance matrix for all the bins used in Figure 4. The covariance matrix from the background-only fit is shown."
table_covariance.location = "Supplementary material"
for var in [x_covariance,y_covariance,z_covariance]:
    table_covariance.add_variable(var)
submission.add_table(table_covariance)
### End covariance


#
outdir = "hig20017_output"
submission.create_files(outdir)



