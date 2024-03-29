
// define a function with 3 parameters
Double_t MyGaus(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) 
        arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

Double_t TotalFit(Double_t *x, Double_t *par){
    return MyGaus(x,par) + MyGaus(x,&par[3]) + MyGaus(x,&par[6]) + MyGaus(x,&par[9]);
}


void PrintStateData(TF1* f,TF1* state,Double_t binWidth){
    Double_t p[3];
    state->GetParameters(&p[0]);

    cout<<p[1]<<", "
        <<f->GetParError(1)<<", "
        <<2.35482*p[2]<<", "
        <<state->Integral(p[1]-4*p[2], p[1]+4*p[2])/binWidth<<", "
        <<sqrt(state->Integral(p[1]-4*p[2], p[1]+4*p[2])/binWidth)<<
        endl;
}


void exc_energy_fit(){

    // Load in rootfile, define and histogram
    // TFile *file_analysis = TFile::Open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_14.root"); 
    // TFile *file_analysis = TFile::Open("/mnt/data0/2014_06_25Al_dn_jbaker/output_2022_10_19.root"); 
    TFile *file_analysis = TFile::Open("C:\\Users\\elite\\OneDrive\\Research\\analysis_25Al-2H-final\\data\\histograms\\output_2023_03_14-1.root"); 

    auto c1 = new TCanvas("c1","",325,325); // 325 pixels is 8.6 cm
    TH1F *h1 = (TH1F*)file_analysis->Get("InvMassExcEnergy");

    gStyle->SetOptStat(0); // turn off statistic box
    h1->SetTitle(";^{26}Si Excitation Energy (MeV); Counts");
    h1->GetXaxis()->SetTitleOffset(1.1);
    //h1->SetFillColorAlpha(kGray, 0.25);
    //h1->SetFillStyle(4050);
    //h1->SetMarkerStyle(20); // Set marker for dots
    //h1->Draw("E1HISTBAR"); // Draw histogram with error bars
    h1->Draw("HIST"); // Draw histogram with error bars
    gStyle->SetErrorX(0); //supresses the x error

    double xmin=5.2,xmax=8.1; //zoom into the graph with these limits
    double ymin=0.0,ymax=45;
    auto axis1 = h1->GetXaxis();
    axis1->SetRangeUser(xmin, xmax); //zoom into graph
    h1->SetAxisRange(ymin,ymax,"Y"); // Set y axis
    gPad->Update(); c1->Update();
    // double ymax=gPad->GetUymax();

    // Make a secondary axis on top for resonance energy
    //xmin = 5.3;
    //xmax = 7.1;
    auto f1=new TF1("f1","x",0,xmax-5.514);
    auto axis2 = new TGaxis(5.514,ymax,xmax,ymax,"f1",510,"S-");
    axis2->SetTitleFont(42);
    axis2->SetLabelFont(42);
    
    axis2->SetTitleSize(0.04);
    axis2->SetLabelSize(0.04);
    
    axis2->SetTickSize(0.04);
    axis2->SetTextFont(42);
    axis2->SetTextSize(0.04);
    
    axis2->SetLabelOffset(-0.01);    

    axis2->SetTitle("Resonance E_{cm} (MeV)");
    axis2->SetLabelFont(gStyle->GetLabelFont("X"));
    axis2->SetLabelSize(gStyle->GetLabelSize("X"));
    axis2->Draw();


    //////////////////////////// Fitting Method //////////////////////////////////////////

    // Define the Global fit
    Int_t NumOfGausPeaks = 4;
    Int_t NumOfParameters = 3*NumOfGausPeaks;
    Double_t param[NumOfParameters]; // there are three parameter for each gaussian

    xmin = 5.6;
    xmax = 8.0;
    auto globalFit = new TF1("globalFit","TotalFit",xmin, xmax, NumOfParameters);

    xmin = 5.6;
    xmax = 6.11;
    auto gaus1 = new TF1("gaus1","MyGaus",xmin,xmax,3); // 3+
    gaus1->SetParLimits(0,10,70);
    gaus1->SetParLimits(1,xmin,xmax);
    gaus1->SetParLimits(2,2E-2,30E-2);

    xmin = 6.0;
    xmax = 6.5;
    auto gaus2 = new TF1("gaus2","MyGaus",xmin,xmax,3); // 3+
    gaus2->SetParLimits(0,10,70);
    gaus2->SetParLimits(1,xmin,xmax);
    gaus2->SetParLimits(2,2E-2,30E-2);

    xmin = 6.5;
    xmax = 7.0;
    auto gaus3 = new TF1("gaus3","MyGaus",xmin,xmax,3); // 3+
    gaus3->SetParLimits(0,10,70);
    gaus3->SetParLimits(1,xmin,xmax);
    gaus3->SetParLimits(2,2E-2,30E-2);

    xmin = 7.0;
    xmax = 8.0;
    auto gaus4 = new TF1("gaus4","MyGaus",xmin,xmax,3); // 3+
    gaus4->SetParLimits(0,10,70);
    gaus4->SetParLimits(1,xmin,xmax);
    gaus4->SetParLimits(2,2E-2,30E-2);

    // Fit all the separted peaks and store their parameters
    h1->Fit(gaus1,"0QR"); // 0=no print, Q=quite no output, R=use user range, + add to fit list
    h1->Fit(gaus2,"0QR+");
    h1->Fit(gaus3,"0QR+");
    h1->Fit(gaus4,"0QR+");


    // Retrieve the local fit parameters and save them into param
    gaus1->GetParameters(&param[0]);
    gaus2->GetParameters(&param[3]);
    gaus3->GetParameters(&param[6]);
    gaus4->GetParameters(&param[9]);
    // Set the global parameters
    globalFit->SetParameters(param);



    // Fit and plot the global function
    // h1->Fit(globalFit,"0MR+");
    auto fitResult = h1->Fit(globalFit,"S0MR+");
    auto covMatrix = fitResult->GetCovarianceMatrix();
    //covMatrix.Print();

    globalFit->GetParameters(&param[0]);
    auto paramError = globalFit->GetParErrors();

    globalFit->SetLineColor(kBlack);
    globalFit->Draw("same C");
    cout<<"---------------------------------------\n";
    cout<<"Reduced Chi-Square = "<<globalFit->GetChisquare()/globalFit->GetNDF()<<endl;
    cout<<"---------------------------------------\n";



    ///////////// Create and plot singular histrograms //////////////////////
    double hBW = h1->GetBinWidth(0);
    double xlow,xhigh;

    gaus1->SetParameters(param[0],param[1],param[2]);
    xlow=param[1]-4.0*param[2];
    xhigh=param[1]+4.0*param[2];

    Int_t col1 = TColor::GetColor("#c8b95d");
    gaus1->SetLineColor(col1);
    gaus1->DrawF1(xlow,xhigh,"same C");

    gaus2->SetParameters(param[3],param[4],param[5]);
    xlow=param[4]-4.0*param[5];
    xhigh=param[4]+4.0*param[5];
    
    Int_t col2 = TColor::GetColor("#2cb4a4");
    gaus2->SetLineColor(col2);
    gaus2->DrawF1(xlow,xhigh,"same C");

    gaus3->SetParameters(param[6],param[7],param[8]);
    xlow=param[7]-4.0*param[8];
    xhigh=param[7]+4.0*param[8];
    
    Int_t col3 = TColor::GetColor("#342c84");
    gaus3->SetLineColor(col3);
    gaus3->DrawF1(xlow,xhigh,"same C");

    gaus4->SetParameters(param[9],param[10],param[11]);
    xlow=param[10]-4.0*param[11];
    xhigh=param[10]+4.0*param[11];

    Int_t col4 = TColor::GetColor("#1484d4");
    gaus4->SetLineColor(col4);
    gaus4->SetLineStyle(2);
    gaus4->DrawF1(xlow,xhigh,"same C");



    ////////////////////  Print Values //////////////////////
    cout<<"peak, peak-error, FWHM, integral, integral-error \n";  
    PrintStateData(globalFit,gaus1,hBW);
    PrintStateData(globalFit,gaus2,hBW);
    PrintStateData(globalFit,gaus3,hBW);
    PrintStateData(globalFit,gaus4,hBW);



    //////////// Add verticals Lines //////////////////
    // Separation energy
    double x0 = 5.514;
    auto l0 = new TLine(x0,0,x0,ymax);  
    l0->SetLineStyle(kDashed); //dashed */
    l0->SetLineColor(kBlack);
    l0->Draw(); 

    //   // double x1 = 5.927; // 3+ */
    //   // auto l1 = new TLine(x1,0,x1,ymax);   
    //   // l1->SetLineStyle(9); //dashed */
    //   // l1->SetLineColor(kGray); //red */
    //   // l1->Draw(); 

    //    double x2 = 6.295; // 2+ 
    //    auto l2 = new TLine(x2,0,x2,ymax);  
    //    l2->SetLineStyle(9); //dashed 
    //    l2->SetLineColor(kGray);  
    //    l2->Draw(); 

    //    double x3 = 6.382; 
    //    auto l3 = new TLine(x3,0,x3,ymax);  
    //    l3->SetLineStyle(9); //dashed 
    //    l3->SetLineColor(kGray);  
    //    l3->Draw(); 

    //    double x4 = 6.46; 
    //    auto l4 = new TLine(x4,0,x4,ymax);  
    //    l4->SetLineStyle(9); //dashed 
    //    l4->SetLineColor(kGray);  
    //    l4->Draw(); 

    //   /* double x5 = 6.787; */
    //   /* auto l5 = new TLine(x5,0,x5,ymax);  */
    //   /* l5->SetLineStyle(9); //dashed */
    //   /* l5->SetLineColor(kGray);  */
    //   /* l5->Draw(); */


    /* auto legend = new TLegend(0.27,0.68,0.85,0.86); */
    /* legend->SetNColumns(2); */
    /* legend->SetBorderSize(0); */
    /* //legend->SetFont(42); */
    /* legend->SetTextSize(0.04); */
    /* // legend->AddEntry(l0,"Sp=5.514","l"); */
    /* legend->AddEntry(globalFit,Form("Global Fit, #chi^{2}_{r} = %.2f",globalFit->GetChisquare()/globalFit->GetNDF()),"l"); */
    /* //legend->AddEntry(gaus1,"Peak 1, 3^{+}_{3} Fit","l"); */
    /* //legend->AddEntry(gaus2,"Peak 2 Fit","l"); //2^{+}_{6} + 2^{+}_{7} */
    /* //legend->AddEntry(gaus3,"Peak 3 Fit","l"); //3^{-}_{1} */
    /* //legend->AddEntry(gaus4,"Background from higher states","l"); */
    /* legend->Draw("SAME"); */







//////////// Add Text //////////////////

//   auto t1 = new TLatex();
//   t1->SetNDC();
//   t1->SetTextFont(42);
//   t1->SetTextSize(0.035);
//   t1->SetTextAngle(90);
//   t1->DrawLatex(0.35,0.63,"5.91 MeV, 3^{+}_{3}");

//   auto t2 = new TLatex();
//   t2->SetNDC();
//   t2->SetTextFont(42);
//   t2->SetTextSize(0.035);
//   t2->SetTextAngle(90);
//   t2->DrawLatex(0.5,0.5,"6.30 MeV, 2^{+}_{6}");

//   auto t3 = new TLatex();
//   t3->SetNDC();
//   t3->SetTextFont(42);
//   t3->SetTextSize(0.035);
//   t3->SetTextAngle(90);
//   t3->DrawLatex(0.70,0.63,"6.56 MeV");

//   auto t4 = new TLatex();
//   t4->SetNDC();
//   t4->SetTextFont(42);
//   t4->SetTextSize(0.035);
//   t4->SetTextAngle(90);
//   t4->DrawLatex(0.75,0.70,"6.72 MeV, 3^{-}_{1}");

//   auto t5 = new TLatex();
//   t4->SetNDC();
//   t4->SetTextFont(42);
//   t4->SetTextSize(0.035);
//   t4->SetTextAngle(90);
//   t4->DrawLatex(0.87,0.4,"7.02 MeV, 3^{+}_{4}");



}
