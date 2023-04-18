#include <stdio.h>
#include <cmath>
#include <vector>
#include "TCanvas.h"

// define functions
template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
void reaction_rate();
void density();
void reaclib();

double res_rate(double t, double wg, double er);
double res_rate_error(double t, double wg, double er, double dwg, double der);
double direct_capture_rate(double t, int zt, double s0);
double reaclib_rate(double t, double coeff[]);

// define global variables
double mu = 24.990428308*1.007825031898/(1.007825031898+24.990428308); // from AME2020

//Temanson et al.
double er_3p = 0.4154; // MeV
double er_3p_error = 0.0008; 
double wg_3p = 2.62E-08; //MeV
double er_wg_3p = wg_3p*0.25;
//Perello et al.
double er_0p = 0.375; //MeV 
double er_0p_error = 0.002;
double wg_0p = 2.20E-10;
double er_wg_0p = wg_0p*0.3;
//Hamil et al.
double er_1p = 0.1622;
double er_1p_error = 0.0003;
double wg_1p = 2.60E-15;
double er_wg_1p = wg_1p*0.3;

double S0 = 0.028;

// temperature
std::vector<double> temp = linspace(0.01,0.71,40);

// the main program, call the three plots
void astro_plots()
{
    reaction_rate();
    density();
    reaclib();
}

// making the plot for the reaction rate
void reaction_rate()
{
    auto c_rate = new TCanvas("c_rate","",325,325);
    c_rate->Draw();
    c_rate->SetLogy();

    std::vector<double> rate_3p;
    std::vector<double> rate_error_3p;
    std::vector<double> rate_0p;
    std::vector<double> rate_error_0p;
    std::vector<double> rate_1p;
    std::vector<double> rate_error_1p;
    std::vector<double> rate_ds;
    std::vector<double> rate_error_ds;

    for(auto& t : temp)
    {
        rate_3p.push_back(res_rate(t,wg_3p,er_3p));
        rate_error_3p.push_back(res_rate_error(t,wg_3p,er_3p, er_wg_3p,er_3p_error));
        rate_0p.push_back(res_rate(t,wg_0p,er_0p));
        rate_error_0p.push_back(res_rate_error(t,wg_0p,er_0p, er_wg_0p,er_0p_error));
        rate_1p.push_back(res_rate(t,wg_1p,er_1p));
        rate_error_1p.push_back(3.0*res_rate_error(t,wg_1p,er_1p, er_wg_1p,er_1p_error));
        rate_ds.push_back(direct_capture_rate(t,13,S0));
        rate_error_ds.push_back(0.3*direct_capture_rate(t,13,S0));
    }

    auto mg = new TMultiGraph();
    auto gr_3p = new TGraphErrors(temp.size(),&temp[0],&rate_3p[0],0,&rate_error_3p[0]);
    auto gr_0p = new TGraphErrors(temp.size(),&temp[0],&rate_0p[0],0,&rate_error_0p[0]);
    auto gr_1p = new TGraphAsymmErrors(temp.size(),&temp[0],&rate_1p[0],0,0,&rate_error_1p[0],0);
    auto gr_dc = new TGraphErrors(temp.size(),&temp[0],&rate_ds[0],0,&rate_error_ds[0]);

    gStyle->SetPalette(kLightTemperature);
    // gStyle->SetPalette(kDarkRainBow);
    mg->Add(gr_3p,"c4");
    mg->Add(gr_0p,"c4");
    mg->Add(gr_dc,"c4");
    mg->Add(gr_1p,"c|>");
    mg->Draw("A plc pfc");
    mg->GetYaxis()->SetNdivisions(506);
    mg->GetYaxis()->SetRangeUser(1E-15,1E1);
    mg->GetXaxis()->SetRangeUser(0,0.6);
    mg->GetYaxis()->SetTitleOffset(1.5);
    mg->GetYaxis()->SetLabelOffset(0.0);
    mg->SetTitle("");
    mg->GetXaxis()->SetTitle("Temperature (GK)");
    mg->GetYaxis()->SetTitle("N_{A}<#sigma#nu> (cm^{3}s^{-1}mol^{-1})");

    auto leg = new TLegend(0.12,0.67,0.48,0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr_3p,"5.9276(10) MeV, 3^{+}_{3}","lf");
    leg->AddEntry(gr_0p,"5.890(8) MeV, 0^{+}_{4}","lf");
    leg->AddEntry(gr_1p,"5.6752(14) MeV, 1^{+}_{1}","lf");
    leg->AddEntry(gr_dc,"direct capture","lf");
    leg->Draw("same");

    // c_rate->Draw("AL");
    c_rate->SaveAs("reaction_rate.png");
    c_rate->SaveAs("reaction_rate.eps");
}

void density()
{
    auto c_density = new TCanvas("c_density","",325,325);
    c_density->Draw();
    c_density->SetLogy();

    std::vector<double> density;
    std::vector<double> density_error;

    double al25_half_life = 7.183/std::log(2.0); // seconds
    double al25_half_life_error = 0.012;

    std::vector<double> temp2 = linspace(0.01,0.71,100);
    for(auto &t : temp2)
    {
        double rate = res_rate(t,wg_3p,er_3p) + res_rate(t,wg_0p,er_0p) + 
            res_rate(t,wg_1p,er_1p) + direct_capture_rate(t,13,S0);
        
        double rate_error = res_rate_error(t,wg_3p,er_3p, er_wg_3p, er_3p_error) + 
            res_rate_error(t,wg_0p,er_0p, er_wg_0p,er_0p_error) + 
            res_rate_error(t,wg_1p,er_1p, er_wg_1p,er_1p_error) + 
            0.3*direct_capture_rate(t,13,S0);

        density.push_back(al25_half_life/rate);
        density_error.push_back(al25_half_life/rate * 
            std::sqrt(
                std::pow(rate_error/rate, 2.0) + 
                std::pow(al25_half_life_error/al25_half_life, 2.0)));
    }

    auto gr_density = new TGraphErrors(temp2.size(),&temp2[0],&density[0]);
    // auto gr_density = new TGraphErrors(temp2.size(),&temp2[0],&density[0],0,&density_error[0]);

    gr_density->GetXaxis()->SetRangeUser(0,0.6);
    gr_density->GetYaxis()->SetRangeUser(1E1,1E15);
    gr_density->GetYaxis()->SetNdivisions(506);
    gr_density->SetTitle(";Temperature (GK);Proton Density (g cm^{-3})");
    gr_density->GetYaxis()->SetTitleOffset(1.35);
    gr_density->GetYaxis()->SetLabelOffset(0.0);
    gr_density->SetLineWidth(2);
    gr_density->Draw("AC");
    
    auto Box = new TBox(0.145,1.0E3,0.418,1.0E4);
    // Box->SetFillColorAlpha(1, 1.5);
    Box->SetFillColor(kBlack);
    gStyle->SetHatchesSpacing(2);
    Box->SetFillStyle(3354);
    Box->Draw("same");

    auto nova = new TLatex();
    nova->DrawLatexNDC(0.35,0.23,"nova");

    auto text1 = new TLatex();
    text1->DrawLatexNDC(0.5, 0.5, "^{25}Al(p,#gamma)^{26}Si");

    auto text2 = new TLatex();
    text2->DrawLatexNDC(0.15, 0.3, "^{25}Al(#beta^{+})");

    // c_density->SaveAs("density.root");
    c_density->SaveAs("density.eps");
    c_density->SaveAs("density.png");
}
//
void reaclib()
{
    auto c_reaclib = new TCanvas("c_reaclib","",325,325);

    double Li2020_coeff1[7] = {-6.207810E+00,-1.731020E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,-2.442940E-01};
    double Li2020_coeff2[7] = {5.387930e+00,-1.245870e+01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.487210e+00};
    double Li2020_coeff3[7] = {8.745920e+00,-4.788620e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.485450e+00};

    double Ma2010_coeff1[7] = {1.981490e+01,0.000000e+00,-2.318660e+01,0.000000e+00,0.000000e+00,0.000000e+00,-6.666670e-01};
    double Ma2010_coeff2[7] = {-7.985670e+00,-1.886890e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    double Ma2010_coeff3[7] = {8.639460e+00,-4.675450e+00,0.000000e+00,-6.012900e-01,2.904360e-01,-4.676690e-03,-1.500000e+00};

    double Il2010_coeff1[7] = {2.311230e+01,0.000000e+00,-2.352420e+01,-9.536160e+00,2.649930e+01,-1.472910e+01,-6.666670e-01};
    double Il2010_coeff2[7] = {-8.132820e+00,-1.896620e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    double Il2010_coeff3[7] = {1.137370e+01,-4.843520e+00,0.000000e+00,-3.389800e+00,1.191190e+00,-1.024290e-01,-1.500000e+00};
    
    std::vector<double> Il2010;
    std::vector<double> Li2020;
    std::vector<double> Ma2010;

    std::vector<double> rate;
    std::vector<double> rate_error;

    for(auto& t : temp)
    {
        double rate_total = 
            res_rate(t,wg_3p,er_3p) + 
            res_rate(t,wg_0p,er_0p) + 
            res_rate(t,wg_1p,er_1p) + 
            direct_capture_rate(t,13,S0);

        Il2010.push_back((reaclib_rate(t,Il2010_coeff1)+reaclib_rate(t,Il2010_coeff2)+reaclib_rate(t,Il2010_coeff3))/rate_total);
        Li2020.push_back((reaclib_rate(t,Li2020_coeff1)+reaclib_rate(t,Li2020_coeff2)+reaclib_rate(t,Li2020_coeff3))/rate_total);
        Ma2010.push_back((reaclib_rate(t,Ma2010_coeff1)+reaclib_rate(t,Ma2010_coeff2)+reaclib_rate(t,Ma2010_coeff3))/rate_total);

        rate.push_back((
            res_rate(t,wg_3p,er_3p) + 
            res_rate(t,wg_0p,er_0p) + 
            res_rate(t,wg_1p,er_1p) + 
            direct_capture_rate(t,13,S0))/rate_total);

        rate_error.push_back((
            res_rate_error(t,wg_3p,er_3p, er_wg_3p, er_3p_error) + 
            res_rate_error(t,wg_0p,er_0p, er_wg_0p, er_0p_error) + 
            res_rate_error(t,wg_1p,er_1p, er_wg_1p, er_1p_error) + 
            0.3*direct_capture_rate(t,13,S0)/rate_total));
    }

    auto gr_Il2010 = new TGraph(temp.size(),&temp[0],&Il2010[0]);
    auto gr_Li2020 = new TGraph(temp.size(),&temp[0],&Li2020[0]);
    auto gr_Ma2010 = new TGraph(temp.size(),&temp[0],&Ma2010[0]);
    auto gr_et2023 = new TGraph(temp.size(),&temp[0],&rate[0]);
    // auto gr_et2023 = new TGraphErrors(temp.size(),&temp[0],&rate[0],0,&rate_error[0]);

    gr_Il2010->SetLineWidth(2);
    gr_Li2020->SetLineWidth(2);
    gr_Ma2010->SetLineWidth(2);
    gr_et2023->SetLineWidth(2);

    auto mg = new TMultiGraph();
    mg->Add(gr_Il2010,"C");
    mg->Add(gr_Li2020,"C");
    mg->Add(gr_Ma2010,"C");
    mg->Add(gr_et2023,"C");
    mg->SetTitle(";Temperature (GK);Ratio of reaction rates");
    mg->GetYaxis()->SetRangeUser(0,6);
    mg->GetXaxis()->SetRangeUser(0,0.6);

    // gStyle->SetPalette(kBird);
    mg->Draw("A plc pfc");

    auto leg = new TLegend(0.15,0.6,0.6,0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr_Il2010,"Iliadis, 2010","l");
    leg->AddEntry(gr_Li2020,"Liang, 2020","l");
    leg->AddEntry(gr_Ma2010,"Matic, 2010","l");
    leg->AddEntry(gr_et2023,"This work","l");
    leg->Draw("same");

    c_reaclib->SaveAs("reaclib.eps");
    c_reaclib->SaveAs("reaclib.png");
}


//  full functions
double res_rate(double t, double wg, double er)
{   
    return 1.5394E11 * wg * std::pow(t*mu,-1.5) * std::exp(-11.605 * er/t);
}
double res_rate_error(double t, double wg, double er, double dwg, double der)
{
    return std::sqrt( 
        std::pow(res_rate(t,wg,er)*dwg/wg, 2.0) + 
        std::pow(res_rate(t,wg,er)*der*(-11.605/t), 2.0) );
}
double direct_capture_rate(double t, int zt, double s0)
{
    return 7.8327E9*std::pow(zt/(std::pow(t,2.0)*mu), 1.0/3.0) * 
        std::exp(-4.2487*std::pow(std::pow(zt,2.0)*mu/t, 1.0/3.0)) * 
        s0*(1+0.09807*std::pow(t/(std::pow(zt,2.0)*mu), 1.0/3.0));
}
double reaclib_rate(double t, double coeff[])
{   
    return std::exp(coeff[0] + 
        coeff[1]/t + 
        coeff[2]/std::pow(t,1.0/3.0) +
        coeff[3]*std::pow(t,1.0/3.0) + 
        coeff[4]*t + 
        coeff[5]*std::pow(t,5.0/3.0) + 
        coeff[6]*std::log(t));
}
template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in)
{
    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1) 
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
    return linspaced;
}
