#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:21:18 2022

@author: est18c
"""

import ROOT as rt

file = rt.TFile.Open("/mnt/data0/2014_06_25Al_dn_jbaker/Filtered_noDS.root","READ")

tree = file.Get("DataTree")

c1 = rt.TCanvas("c1")

tree.Draw("ADC3[0]")

c1.Update()