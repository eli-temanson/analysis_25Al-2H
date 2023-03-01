#!/usr/bin/env python
# coding: utf-8

# # PyROOT plotting for the $^{25}\textrm{Al}(d,n)^{26}\textrm{Si}$ paper

import ROOT
import numpy as np

col1 = ROOT.TColor.GetColor("#c8b95d")
col2 = ROOT.TColor.GetColor("#2cb4a4")
col3 = ROOT.TColor.GetColor("#342c84")
col4 = ROOT.TColor.GetColor("#1484d4")

mu=24.990428308*1.007825031898/(1.007825031898+24.990428308) # from AME2020
temp = np.linspace(0.06,1.0,100)
temp_short = np.linspace(0.06,1.0,35)

def Rate(t, wg, er):
    return (1.5394E11)*wg*(t*mu)**(-1.5)*np.exp(-11.605*er/t)

def RateError(t,wg,er,dwg,der):
    return np.sqrt( (Rate(t,wg,er)*dwg/wg)**2 + (Rate(t,wg,er)*der*(-11.605/t))**2 ) 

def DCRate(t,zt,s0):
    return (7.8327E9)*(zt/((t**2)*mu))**(1/3)*np.exp(-4.2487*((zt**2)*mu/t)**(1/3))*s0*(1+0.09807*(t/((zt**2)*mu))**(1/3))

# ## Plot the reaction rate for both reasonance and direct capture

Er_3p=0.4154 # MeV
Er_3p_error=0.0008
wg_3p=1.78E-08 # MeV

Er_0p=0.3671
Er_0p_error=0.0003
wg_0p=2.40E-10

Er_1p=0.1622
Er_1p_error=0.0003
wg_1p=2.60E-15

S0 = 0.028

c1 = ROOT.TCanvas('c1', '', 325, 325)
c1.Draw()

gr_3p=ROOT.TGraphErrors(temp.size,temp,Rate(temp,wg_3p,Er_3p), 0,RateError(temp,wg_3p,Er_3p, wg_3p*0.25,Er_3p_error))
gr_0p=ROOT.TGraphErrors(temp.size,temp,Rate(temp,wg_0p,Er_0p), 0,RateError(temp,wg_0p,Er_0p, wg_0p*0.3,Er_0p_error))
# gr_1p=ROOT.TGraphErrors(temp.size,temp,Rate(temp,2.60E-15,0.1614), 0,RateError(temp,2.60E-15,0.1614, 2.60E-15*0.3,0.0015))
gr_1p=ROOT.TGraphAsymmErrors(temp_short.size,temp_short,Rate(temp_short,wg_1p,Er_1p), 0,0, 2.5*RateError(temp_short,wg_1p,Er_1p, wg_1p*0.3,Er_1p_error),0)

gr_dc=ROOT.TGraphErrors(temp.size,temp,DCRate(temp,13,S0), 0, DCRate(temp,13,S0)*0.3)

Total_Rate = Rate(temp,wg_3p,Er_3p)+Rate(temp,wg_0p,Er_0p)+Rate(temp,wg_1p,Er_1p)+DCRate(temp,13,S0)
Total_Rate_Error = np.sqrt(RateError(temp,wg_3p,Er_3p, wg_3p*0.25,Er_3p_error)**2 +  RateError(temp,wg_0p,Er_0p, wg_0p*0.3,Er_0p_error)**2 + RateError(temp,wg_1p,Er_1p, wg_1p*0.3,Er_1p_error)**2 + (DCRate(temp,13,S0)*0.3)**2) 

gr_total_rate=ROOT.TGraphErrors(temp.size,temp,Total_Rate, 0,Total_Rate_Error)

# ROOT.gStyle.SetPalette(57,0,0.7) # kBird

gr_total_rate.SetTitle('')
gr_total_rate.GetXaxis().SetTitle('Temperature (Gk)')
gr_total_rate.GetYaxis().SetTitle('N_{A}<#sigma#nu> (cm^{3}s^{-1}mol^{-1})')
gr_total_rate.GetXaxis().SetRangeUser(0.09,0.7)
gr_total_rate.GetYaxis().SetNdivisions(508)
gr_total_rate.GetYaxis().SetTitleOffset(1.65)

# gr_total_rate.SetLineWidth(2)
gr_total_rate.SetLineStyle(1)
gr_total_rate.SetLineColor(1)
gr_total_rate.SetFillColor(1)

# gr_3p.SetLineWidth(2)
# gr_3p.SetLineStyle(1)
# gr_0p.SetLineWidth(2)
# gr_0p.SetLineStyle(2)
# gr_1p.SetLineWidth(2)
# gr_1p.SetLineStyle(1)
# gr_dc.SetLineWidth(2)
# gr_dc.SetLineStyle(4)

gr_3p.SetFillStyle(1001)
gr_0p.SetFillStyle(1001)
gr_dc.SetFillStyle(1001)
gr_1p.SetFillStyle(1001)

gr_3p.SetLineColor(col1)
gr_0p.SetLineColor(col2)
gr_dc.SetLineColor(col3)
gr_1p.SetLineColor(col4)

gr_3p.SetFillColor(col1)
gr_0p.SetFillColor(col2)
gr_dc.SetFillColor(col3)
gr_1p.SetFillColor(col4)

gr_total_rate.Draw('A4 same')
gr_3p.Draw('A4 same')
gr_0p.Draw('A4 same')
gr_dc.Draw('A4 same')
gr_1p.Draw('CL> same')

c1.SetLogy()

leg = ROOT.TLegend(0.15,0.7,0.48,0.9)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
# leg.SetTextFont(42)
# leg.SetTextSize(0.035)
leg.AddEntry(gr_total_rate,'Total Rate','lf')
leg.AddEntry(gr_3p,'5.9276(10) MeV, 3^{+}_{3}','lf')
leg.AddEntry(gr_0p,'5.890(8) MeV, 0^{+}_{4}','lf')
leg.AddEntry(gr_1p,'5.6752(14) MeV, 1^{+}_{1}','lf')
leg.AddEntry(gr_dc,'Direct capture','lf')
leg.Draw("same")

c1.SetLeftMargin(0.13)
c1.SaveAs("rrate.png")
# c1.SaveAs("rrate.root")

# ==========================================================

# Plot the density-temperature profile where the (p,γ) lifetime dominates the β+

Al25_Half_Life=7.183/np.log(2) #seconds
Al25_Half_Life_Error=0.012

Density_Temp_Error = (Al25_Half_Life/Total_Rate)*np.sqrt((Total_Rate_Error/Total_Rate)**2 + (Al25_Half_Life_Error/Al25_Half_Life)**2 )

c2 = ROOT.TCanvas('c2', '', 325, 325)
Density_Temp=ROOT.TGraphErrors(temp.size,temp,Al25_Half_Life/Total_Rate, 0, Density_Temp_Error)

c2.Draw()
ROOT.gStyle.SetPalette(52,0,0.7) # kgreyScale

Density_Temp.Draw('AC4 PLC PFC')
Density_Temp.SetTitle('')
Density_Temp.GetXaxis().SetTitle('Temperature (Gk)')
Density_Temp.GetYaxis().SetTitle('Proton Density (g cm^{-3})')
Density_Temp.GetXaxis().SetRangeUser(0,0.7)
Density_Temp.GetYaxis().SetNdivisions(508)
Density_Temp.GetYaxis().SetTitleOffset(1.45)

Box=ROOT.TBox(0.145,1.0E3,0.418,1.0E4)
Box.SetFillColorAlpha(1, 1.5)
Box.SetFillStyle(3004)
Box.Draw('same')

c2.SetLogy()

text1=ROOT.TLatex()
text1.DrawLatexNDC(.5,.5,"^{25}Al(p,#gamma)^{26}Si")

text2=ROOT.TLatex()
text2.DrawLatexNDC(.25,.25,"^{25}Al(#beta^{+})")

# text3=ROOT.TLatex();
# text3.DrawLatexNDC(.25,.75,"#it{#tau_{p} = #tau_{#beta}}")

text4=ROOT.TLatex()
text4.DrawLatexNDC(.50,.20,"nova")

c2.SetLeftMargin(0.13)
c2.SaveAs("density.png")
# c2.SaveAs("density.root")

# ==========================================================

# ## Plotting how the reaction rate compares to that of the suggest JINA REACLIB database

def ReaclibRate(t,a):
    return np.exp(a[0] + a[1]/t + a[2]/(t**(1/3))+ a[3]*t**(1/3) + a[4]*t + a[5]*t**(5/3) + a[6]*np.log(t))

Li2020_coeff1 = [-6.207810E+00,-1.731020E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,-2.442940E-01]
Li2020_coeff2 = [5.387930e+00,-1.245870e+01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.487210e+00]
Li2020_coeff3 = [8.745920e+00,-4.788620e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.485450e+00]
Li2020_rate=ReaclibRate(temp,Li2020_coeff1)+ReaclibRate(temp,Li2020_coeff2)+ReaclibRate(temp,Li2020_coeff3)

Ma2010_coeff1 = [1.981490e+01,0.000000e+00,-2.318660e+01,0.000000e+00,0.000000e+00,0.000000e+00,-6.666670e-01]
Ma2010_coeff2 = [-7.985670e+00,-1.886890e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00]
Ma2010_coeff3 = [8.639460e+00,-4.675450e+00,0.000000e+00,-6.012900e-01,2.904360e-01,-4.676690e-03,-1.500000e+00]
Ma2010_rate=ReaclibRate(temp,Ma2010_coeff1)+ReaclibRate(temp,Ma2010_coeff2)+ReaclibRate(temp,Ma2010_coeff3)

Il2010_coeff1 = [2.311230e+01,0.000000e+00,-2.352420e+01,-9.536160e+00,2.649930e+01,-1.472910e+01,-6.666670e-01]
Il2010_coeff2 = [-8.132820e+00,-1.896620e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00]
Il2010_coeff3 = [1.137370e+01,-4.843520e+00,0.000000e+00,-3.389800e+00,1.191190e+00,-1.024290e-01,-1.500000e+00]
Il2010_rate=ReaclibRate(temp,Il2010_coeff1)+ReaclibRate(temp,Il2010_coeff2)+ReaclibRate(temp,Il2010_coeff3)

c3 = ROOT.TCanvas('c3', '', 325, 325)
c3.Draw()

gr_li2020=ROOT.TGraph(temp.size,temp,100*(Li2020_rate-Total_Rate)/Total_Rate)
gr_ma2010=ROOT.TGraph(temp.size,temp,100*(Ma2010_rate-Total_Rate)/Total_Rate)
gr_il2010=ROOT.TGraph(temp.size,temp,100*(Il2010_rate-Total_Rate)/Total_Rate)
gr_present=ROOT.TGraphErrors(temp.size,temp,100*(Total_Rate-Total_Rate)/Total_Rate, 0,100*Total_Rate_Error/Total_Rate)

# ROOT.gStyle.SetPalette(57,0,0.7) # kBird
# ROOT.TColor.InvertPalette()
col1 = ROOT.TColor.GetColor("#c8b95d")
col2 = ROOT.TColor.GetColor("#2cb4a4")
col3 = ROOT.TColor.GetColor("#342c84")
col4 = ROOT.TColor.GetColor("#1484d4")

gr_li2020.SetLineColor(col1)
gr_ma2010.SetLineColor(col2)
gr_il2010.SetLineColor(col3)
gr_present.SetFillColor(col4)

gr_il2010.SetTitle('')
gr_il2010.GetXaxis().SetTitle('Temperature (Gk)')
gr_il2010.GetYaxis().SetTitle('Reaction Rate Percent Difference (%)')
gr_il2010.GetYaxis().SetTitleOffset(1.65)
gr_il2010.GetXaxis().SetRangeUser(0.09,0.7)
gr_il2010.GetYaxis().SetRangeUser(-60,350)

# gr_li2020.SetLineWidth(4)
# gr_ma2010.SetLineWidth(4)
# gr_il2010.SetLineWidth(4)

gr_il2010.Draw('ACL LC')
gr_li2020.Draw(' CL LC')
gr_ma2010.Draw(' CL LC')

gr_present.SetFillStyle(3325)
gr_present.Draw('CL4 LC FC') #CL4 PLC PFC

# text5=ROOT.TLatex();
# text5.DrawLatexNDC(.50,.20,'')

leg2 = ROOT.TLegend(0.15,0.7,0.48,0.9)
leg2.SetBorderSize(0)
leg2.SetFillColor(0)
leg2.SetFillStyle(0)
leg2.SetTextFont(42)
leg2.SetTextSize(0.035)
leg2.SetNColumns(2)
leg2.AddEntry(gr_il2010,'Iliadis, 2010','l')
leg2.AddEntry(gr_li2020,'Liang, 2020','l')
leg2.AddEntry(gr_ma2010,'Matic, 2010','l')
leg2.AddEntry(gr_present,'Present Work','lf')
leg2.Draw()

c3.SetLeftMargin(0.13)
# c3.SaveAs("difference.root")
c3.SaveAs("difference.png")
