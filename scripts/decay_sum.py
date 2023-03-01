
# python -i ../../py-scripts/decay_sum.py

import ROOT as rt
rfile = rt.TFile.Open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_10.root","READ")
rtree = rfile.Get("ProcessTree")
c1 = rt.TCanvas("c1")
rtree.Draw("decay_sum >> TH1F(48,79,116)","ex_energy < 6.15 && ex_energy > 5.65","")

# 

import uproot as up
import pandas  as pd
import matplotlib.pyplot as plt

utree = up.open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_10.root:ProcessTree")
data = utree.arrays(["decay_sum"],"(ex_energy < 6.15) & (ex_energy > 5.65)", library="pd")

print(data.head())

hist = data.hist(bins=48,range=[79,116])
plt.show()

data.to_csv('~/25Al+d_2014/R-scripts/decay_sum_3p.csv')  