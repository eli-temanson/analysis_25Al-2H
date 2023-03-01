
# import ROOT as rt
# rfile = rt.TFile.Open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_06.root","READ")
# rtree = rfile.Get("ProcessTree")
# c1 = rt.TCanvas("c1")
# rtree.Draw("ex_energy >> TH1F(120,4.514,9.514)")

import uproot as up
import pandas  as pd
import matplotlib.pyplot as plt

utree = up.open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_06.root:ProcessTree")
data = utree.arrays(["ex_energy"], library="pd")

print(data.head())

# plt.figure()
hist = data.hist(bins=120,range=[4.514, 9.514])
plt.show()

data.to_csv('../R-scripts/ex_energy.csv')  