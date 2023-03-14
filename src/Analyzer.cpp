
#ifndef __Analyzer__cpp
#define __Analyzer__cpp

#include "Analyzer.hpp"

Analyzer::Analyzer(TChain *chain,TFile* file,THashTable *hash){
    data_chain = chain;
    out_file = file;
    table0 = hash;
    table0->SetOwner(false);

    // Plotting and running flags
    saveTree=false;
    plot_si_flag=true;
    plot_ic_flag=true;
    plot_timing_flag=true;
    plot_physics_flag=true;
    plot_run_flag=false;
    plot_neutron_flag=false;
    plot_summary_flag=true;

    std::cout<<"Initalizing Detectors \n";
    S2.init("S2",16,16);
    S1.init("S1",16,16);
    ic.init("ic");
    neut.init("neut",12);
    rf.init("rf");
    MCP.init("MCP");
    IonCds.init("IonCds");
    S1Trig.init("S1Trig");

    // This is used in down-scaled reconstruction!
    target.density(0.7987); // g/cm3
    target.thickness(0.000516); // g/cm2
    kapton.density(1.42);
    kapton.thickness(0.000994);
    ic_dl.density(0.00011);
    ic_dl.thickness(0.00036);
    ion_chamber.add(ic_dl);
    ic_x.density(0.00011);
    ic_x.thickness(0.00044);
    ion_chamber.add(ic_x);
    ic_y.density(0.00011);
    ic_y.thickness(0.00044);
    ion_chamber.add(ic_y);
}


void Analyzer::SetReactionInfo(){ std::cout<<"Setting Reaction Info \n";
    rxn1.SetReaction("25Al","2H","n","26Si","1H","25Al");
    rxn1.SetICReco(0.0,0.843,30.2);
    rxn1.SetBeamKE(98.9); 
    rxn1.SetBeamSigma(0.212);
    rxn1.SetFragEst(95.862); // estimate fragment TKE at 0-degrees (from kinematics for choosen excitation energy)
}


void Analyzer::LoadInputData(){ std::cout<<"Loading Input Data \n";
    VarMap.LoadParams("assets/input/Silicon.in");
    VarMap.LoadParams("assets/input/SiliconTime.in");
    VarMap.LoadParams("assets/input/S1FrontGMatch.in");
    VarMap.LoadParams("assets/input/S2FrontGMatch.in");
    VarMap.LoadParams("assets/input/S1BackGMatch.in");
    VarMap.LoadParams("assets/input/S2BackGMatch.in");
    VarMap.LoadParams("assets/input/SiliconChannelMap.in"); 
    VarMap.LoadParams("assets/input/IonChamber.in");
    VarMap.LoadParams("assets/input/IonChamberGridGain.in");
    VarMap.LoadParams("assets/input/Neutron.in");
    VarMap.LoadParams("assets/input/Timing.in"); 
}

void Analyzer::LoadGates(){ std::cout<<"Reading/Loading Gates \n";
    ReadGates("assets/gates/si_protons_small.root");
    ReadGates("assets/gates/si_theta_corr_tight.root");
    ReadGates("assets/gates/si_phi_corr.root");
    ReadGates("assets/gates/ic_band1_2022.root");
    ReadGates("assets/gates/ds_al25_band.root");
    ReadGates("assets/gates/ic_cluster_corr_tight.root");
    ReadGates("assets/gates/ic_t_Et.root");
    ReadGates("assets/gates/ic_p1.root");
    ReadGates("assets/gates/ic_p2.root");
    ReadGates("assets/gates/ic_p3.root");
    ReadGates("assets/gates/InvMassExcEnergy_vFrag.root");


}

void Analyzer::ReadGates(std::string input){
    TFile in(input.c_str());   
    if(in.Get("protons")){ in.GetObject("protons",protons); }
    if(in.Get("ds_al25_band")){ in.GetObject("ds_al25_band",ds_al25_band); }
    if(in.Get("si_theta_corr")){ in.GetObject("si_theta_corr",si_theta_corr); }
    if(in.Get("si_phi_corr")){ in.GetObject("si_phi_corr",si_phi_corr); }
    if(in.Get("ic_band1")){ in.GetObject("ic_band1",ic_band1); }
    if(in.Get("ic_band2")){ in.GetObject("ic_band2",ic_band2); }
    if(in.Get("ic_cluster_corr")){ in.GetObject("ic_cluster_corr",ic_cluster_corr); }
    if(in.Get("ic_t_Et")){ in.GetObject("ic_t_Et",ic_t_Et); }
    if(in.Get("ic_p1")){ in.GetObject("ic_p1",ic_p1); }
    if(in.Get("ic_p2")){ in.GetObject("ic_p2",ic_p2); }
    if(in.Get("ic_p3")){ in.GetObject("ic_p3",ic_p3); }
    if(in.Get("InvMassExcEnergy_vFrag")){ in.GetObject("InvMassExcEnergy_vFrag",InvMassExcEnergy_vFrag);}
    in.Close();
}

void Analyzer::SetBranches(){ std::cout<<"Setting Root Branches \n";
    data_chain->SetBranchAddress("ADC1", &adc1); //Pointers to raw branches
    data_chain->SetBranchAddress("ADC2", &adc2);
    data_chain->SetBranchAddress("ADC3", &adc3);
    data_chain->SetBranchAddress("ADC4", &adc4);
    data_chain->SetBranchAddress("ADC5", &adc5);
    data_chain->SetBranchAddress("ADC6", &adc6);
    data_chain->SetBranchAddress("TDC1", &tdc1);
    data_chain->SetBranchAddress("TDC2", &tdc2);
    data_chain->SetBranchAddress("TDC3", &tdc3);
    data_chain->SetBranchAddress("QDC3", &qdc3);
    if(plot_run_flag){
        data_chain->SetBranchAddress("RunNum", &runNum);
    }
    // Declare all your Branches you want!
    // Save all the raw channel information to be loaded for reevaluation if needed
    // It can also be used for filtering data so the loop is much shorter
    // DataTree->Branch("ADC1",&adc1,"ADC1[32]/F");
    // DataTree->Branch("ADC2",&adc2,"ADC2[32]/F");
    // DataTree->Branch("ADC3",&adc3,"ADC3[32]/F");
    // DataTree->Branch("ADC4",&adc4,"ADC4[32]/F");
    // DataTree->Branch("ADC5",&adc5,"ADC5[32]/F");
    // DataTree->Branch("ADC6",&adc6,"ADC6[32]/F");
    // DataTree->Branch("TDC1",&tdc1,"TDC1[32]/F");
    // DataTree->Branch("TDC2",&tdc2,"TDC2[32]/F");
    // DataTree->Branch("TDC3",&tdc3,"TDC3[32]/F");
    // DataTree->Branch("QDC3",&qdc3,"QDC3[32]/F");
    // if(plot_run_flag){
    //   DataTree->Branch("RunNum",&runNum,"RunNum/I");
    // }

    //CleanTree->Branch("ic_E",&adc4[13],"ic_E/F");
    //CleanTree->Branch("ic_dE",&adc4[12],"ic_dE/F");
    //ProcessTree->Branch("ex_energy",&inv_mass_ex, "ex_energy/F");
    //ProcessTree->Branch("decay_sum",&decay_sum, "decay_sum/F");

}


void Analyzer::SetCalibrations(){ std::cout<<"Setting Detector Calibrations \n";
    ic.SetCalibrations(VarMap); neut.SetCalibrations(VarMap);
    S1.SetCalibrations(VarMap); S2.SetCalibrations(VarMap);
    rf.SetCalibrations(VarMap); 
}


// Insert Your raw data to detector class
void Analyzer::GetEventEntry(){
    S2.Reset(); S1.Reset(); ic.Reset(); neut.Reset();
    double ch;
    for(int i = 0; i<16; i++)
    {
        ch = static_cast<double>(i);
        //Detector 1, Si-dE
        S1.Front.InsertHit(adc2[i+16],0.0,ch);//ring
        S1.Back.InsertHit(adc2[i],tdc2[i],ch);//segment/wedge
        //Detector 2, Si-E
        S2.Front.InsertHit(adc3[i+16],0,ch);//ring
        S2.Back.InsertHit(adc3[i],tdc2[i+16],ch);//segment/wedge  

        //There are 32 channels for each x/y grid
        ic.xgrid.InsertHit(adc5[i],0,ch);
        ic.xgrid.InsertHit(adc6[i+16],0.,ch+16.0);
        ic.ygrid.InsertHit(adc6[i],0,ch);
        ic.ygrid.InsertHit(adc5[i+16],0.,ch+16.0);
        //Insert neutron information
        // if(i < 12){
        //   neut.det.InsertHit(qdc3[i+16],qdc3[i],tdc2[i+16],i);
        // }
    }
    //Insert single channel detector information
    ic.ESeg.InsertHit(adc4[13],tdc1[1],0.);
    ic.DESeg.InsertHit(adc4[12],0.,0.);

    rf.InsertHit(tdc1[0]);
    MCP.InsertHit(tdc1[2]);
    S1Trig.InsertHit(tdc1[3]); //trigger 
    IonCds.InsertHit(tdc1[4]); //trigger downscaled by 1000 or 2000? 

    S1.Front.SortByEnergy(); S1.Back.SortByEnergy();
    S2.Front.SortByEnergy(); S2.Back.SortByEnergy();
    ic.xgrid.SortByEnergy(); ic.ygrid.SortByEnergy();
    // neut.det.SortByEnergy();

} // END OF GET EVENT ENTRY


// =========================================================


void Analyzer::ProcessEvent(){
    // If you didn't create run number data object in evt2root converter you can do it here
    // std::string fname = data_chain->GetFile()->GetName();
    // std::size_t loc = fname.find("run");
    // runNum=std::stoi(fname.substr(loc+3,loc+7)); //string to int, taking 4 digits

    //////////////////////////////////////////////
    // Monitoring events are done here
    if(false)
    // if(true)
    {
        if(!InsideGate(IonCds.TRaw(),2500,2510)) return; 
        if(!InsideGate(ic.ESeg.TRaw(),391,404)) return; 
        if(InsideGate(S1.Back.TRaw(),1,4096)) return; 
        if(InsideGate(S2.Back.TRaw(),1,4096)) return; 

        MyFill(table0,"DS_ic_xGrid_mult",32,0,32,ic.xgrid.Mult());
        MyFill(table0,"DS_ic_yGrid_mult",32,0,32,ic.ygrid.Mult());

        ic.ReconstructCluster(ic.xgrid,ic.xcluster);
        ic.ReconstructCluster(ic.ygrid,ic.ycluster);

        // pile-up rejections
        if(ic.xcluster.Mult()>=2 || ic.ycluster.Mult()>=2) return;
        if(!InsideGate(ic.ycluster.E()-ic.xcluster.E(),591,1450)) return; 
        if(!InsideGate(ic.xcluster.E(),886,1482)) return; 
        
        MyFill(table0,"DS_ic_xClust_mult",32,0,32,ic.xcluster.Mult());
        MyFill(table0,"DS_ic_yClust_mult",32,0,32,ic.ycluster.Mult());
        MyFill(table0,"DS_ic_xgrid_ERaw",4096,0,4096,ic.xgrid.ERaw());
        MyFill(table0,"DS_ic_ygrid_ERaw",4096,0,4096,ic.ygrid.ERaw());
        MyFill(table0,"DS_ic_xGridE_yGridE",512,0,4096,ic.xgrid.E(),512,0,4096,ic.ygrid.E());
        MyFill(table0,"DS_ic_xGridE_dESeg",512,0,4096,ic.xgrid.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"DS_ic_xGridE_ESeg",512,0,4096,ic.xgrid.E(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"DS_ic_yGridE_dESeg",512,0,4096,ic.ygrid.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"DS_ic_yGridE_ESeg",512,0,4096,ic.ygrid.E(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"DS_ic_xClustE",1024,0,4096,ic.xcluster.E());  
        MyFill(table0,"DS_ic_yClustE",1024,0,4096,ic.ycluster.E());
        MyFill(table0,"DS_ic_xClustE_yClustE",512,0,4096,ic.xcluster.E(),512,0,4096,ic.ycluster.E());
        
        MyFill(table0,"DS_ic_xClustE-yClustE_xClustE",512,0,4096,ic.xcluster.E(),1024,-4096,4096,ic.ycluster.E()-ic.xcluster.E());
        
        MyFill(table0,"DS_ic_xClustE_dESeg",512,0,4096,ic.xcluster.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"DS_ic_xClustE_ESeg",512,0,4096,ic.xcluster.E(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"DS_ic_yClustE_dESeg",512,0,4096,ic.ycluster.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"DS_ic_yClustE_ESeg",512,0,4096,ic.ycluster.E(),1000,0,100,ic.ESeg.E());

        MyFill(table0,"DS_IonC_RawPID",2048,0,2048,ic.ESeg.ERaw(),2048,0,2048,ic.DESeg.ERaw());    
        MyFill(table0,"DS_IonC_PID_0",10000,0,100,ic.ESeg.E(),5000,0,50,ic.DESeg.E());
        MyFill(table0,"DS_IonC_PID_1",10000,0,100,ic.ESeg.E()+ic.DESeg.E(),5000,0,50,ic.DESeg.E());
        MyFill(table0,"DS_IonC_ESeg_ERaw",4096,0,4096,ic.ESeg.ERaw());
        MyFill(table0,"DS_IonC_DESeg_ERaw",4096,0,4096,ic.DESeg.ERaw());

        if(InsideGate(ic.ESeg.T(),118,121)){
            MyFill(table0,"DS_IonC_RawPID_t1",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());    
            MyFill(table0,"DS_IonC_PID_t1",6000,0,120,ic.ESeg.E()+ic.DESeg.E(),2500,0,50,ic.DESeg.E());
        }
        if(InsideGate(ic.ESeg.T(),121,124)){
            MyFill(table0,"DS_IonC_RawPID_t2",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());    
            MyFill(table0,"DS_IonC_PID_t2",6000,0,120,ic.ESeg.E()+ic.DESeg.E(),2500,0,50,ic.DESeg.E());
        }
        if(InsideGate(ic_p1,ic.ESeg.ERaw(),ic.DESeg.ERaw())) // selection the 25Al beam
        {
            double ebeam = ic.ESeg.E() + ic.DESeg.E();
            ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_y);
            ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_x);
            ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_dl);
            ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),kapton);
            ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target);
            MyFill(table0,"DS_IonC_beam_reco_p1",600,0,120,ebeam);

            MyFill(table0,"DS_IonC_ESeg_Raw_p1",4096,0,4096,ic.ESeg.ERaw());
            MyFill(table0,"DS_IonC_DESeg_Raw_p1",4096,0,4096,ic.DESeg.ERaw());
            MyFill(table0,"DS_IonC_E_p1",200,0,100,ic.ESeg.E());
            MyFill(table0,"DS_IonC_dE_p1",200,0,100,ic.DESeg.E());
            MyFill(table0,"DS_IonC_RawPID_p1",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());    
        }
        if(InsideGate(ic_p2,ic.ESeg.ERaw(),ic.DESeg.ERaw()))
        {
            MyFill(table0,"DS_IonC_ESeg_Raw_p2",4096,0,4096,ic.ESeg.ERaw());
            MyFill(table0,"DS_IonC_DESeg_Raw_p2",4096,0,4096,ic.DESeg.ERaw());
            MyFill(table0,"DS_IonC_E_p2",200,0,100,ic.ESeg.E());
            MyFill(table0,"DS_IonC_dE_p2",200,0,100,ic.DESeg.E());
            MyFill(table0,"DS_IonC_RawPID_p2",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());   
        }
        if(InsideGate(ic_p3,ic.ESeg.ERaw(),ic.DESeg.ERaw()))
        {
            MyFill(table0,"DS_IonC_ESeg_Raw_p3",4096,0,4096,ic.ESeg.ERaw());
            MyFill(table0,"DS_IonC_DESeg_Raw_p3",4096,0,4096,ic.DESeg.ERaw());
            MyFill(table0,"DS_IonC_E_p3",200,0,100,ic.ESeg.E());
            MyFill(table0,"DS_IonC_dE_p3",200,0,100,ic.DESeg.E());
            MyFill(table0,"DS_IonC_RawPID_p3",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());   
        }

        MyFill(table0,"DS_IonCt",8192,-1024,1024,ic.ESeg.T());
        MyFill(table0,"DS_IonC_rf",8192,0,1024,IonCds.T()-rf.T());
        MyFill(table0,"DS_IonC_rf_mod",8192,-1024,1024,IonCds.T()-rf.TMod());

        if(plot_run_flag)
        {
            MyFill(table0,"DS_run_ic_ESeg_ERaw",200,3200,3400,runNum,4096,0,4096,ic.ESeg.ERaw());
            MyFill(table0,"DS_run_ic_DESeg_ERaw",200,3200,3400,runNum,4096,0,4096,ic.DESeg.ERaw());
            MyFill(table0,"DS_run_ic_ESeg_TRaw",200,3200,3400,runNum,4096,0,4096,ic.ESeg.TRaw());
            MyFill(table0,"DS_run_ic_E",200,3200,3400,runNum,200,0,100,ic.ESeg.E());
            MyFill(table0,"DS_run_ic_DE",200,3200,3400,runNum,200,0,100,ic.DESeg.E());
        }

        // if(InsideGate(ds_al25_band, ic.ESeg.ERaw(), ic.DESeg.ERaw()))
        // {
        //   MyFill(table0,"DS_Al25_ic_E",200,0,100,ic.ESeg.E());
        //   MyFill(table0,"DS_Al25_ic_DE",200,0,100,ic.DESeg.E());
        //   MyFill(table0,"DS_IonC_RawAl25",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());
        //   MyFill(table0,"DS_IonC_E+dE_Al25",200,0,200,ic.ESeg.E() + ic.DESeg.E());
        //   MyFill(table0,"DS_IonC_Theta",100,0,5,ic.Theta()*TMath::RadToDeg());

        //   MyFill(table0,"DS_Al25_beam_reco",70,75,110,efinal);
        //   MyFill(table0,"DS_Al25_beam_reco",60,95,110,ereco);
        //   if(saveTree) DataTree->Fill(); /*Fill the TTree*/ 
        // }
        savedEvents+=1;
    }

    // if(false)
    //     // Check anything inside silicon trigger OR anything outside IC down-scaled trigger
    //     // if(!InsideGate(IonCds.TRaw(),2490,2510) || InsideGate(S1Trig.TRaw(),1,4096)) 
    // {
    //     // MyFill(table0,"Si_PID_Back_nG",1200,0,60,S2.Back.E()+S1.Back.E(), 600,0,20,S2.Back.E());
    //     // MyFill(table0,"Si_PID_Front_nG",1200,0,60,S2.Front.E()+S1.Front.E(), 600,0,20,S2.Front.E());

    //     if(!InsideGate(S1Trig.T()-rf.TMod(), 547., 582.)) return;
    //     if(!InsideGate(ic.ESeg.T(), 200, 600)) return; // 800
    //     // if(!InsideGate(ic_t_Et, ic.ESeg.T(), ic.ESeg.E()+ic.DESeg.E())) return; // no longer needed (ET, 2022)
    //     if(!InsideGate(ic_band1, ic.ESeg.ERaw(), ic.DESeg.ERaw())) return; 
    //     if(!InsideGate(protons, S2.Back.E()+S1.Back.E(), S2.Back.E())) return;

    //     ic.ReconstructCluster(ic.xgrid,ic.xcluster);
    //     ic.ReconstructCluster(ic.ygrid,ic.ycluster);
    //     if(ic.xcluster.Mult()>=2) return;
    //     if(ic.ycluster.Mult()>=2) return;
    //     // if(!InsideGate(ic_cluster_corr,ic.xcluster.E(),ic.ycluster.E())) return; // no longer needed (ET, 2021)

    //     // ic.ReconstructHitPos(ic.xgrid.Ch(),ic.ygrid.Ch()); // obselete, updated to using clusters (ET, 2021) 
    //     ic.ReconstructHitPos(ic.xcluster.Ch(),ic.ycluster.Ch());
    //     S2.ReconstructClusters(S2); S1.ReconstructClusters(S1); 

    //     // Gate on stuff after the reconstruction!
    //     // if(!InsideGate(protons,S2.E()+S1.E(),S2.E())) return;
    //     if(!InsideGate(si_theta_corr, S2.Theta()*TMath::RadToDeg(), S1.Theta()*TMath::RadToDeg())) return; 
    //     if(!InsideGate(si_phi_corr, S2.Phi()*TMath::RadToDeg(), S1.Phi()*TMath::RadToDeg())) return;    

    //     // Calculate rxn1 stuff and then fill histograms
    //     // if(plot_physics_flag){
    //     rxn1.DoPhysics(S1,S2,ic);
    //     // }
    //     // if(rxn1.Fragment.KE < 79.0) return; 
    //     // if(!InsideGate(InvMassExcEnergy_vFrag, rxn1.InvMassExcEnergy,rxn1.Fragment.KE)) return;
    //     inv_mass_ex=rxn1.InvMassExcEnergy;
    //     decay_sum=rxn1.Decay_Heavy.KE+rxn1.Decay_Light.KE;
    //     // FinalTree->Fill();

    //     // if(InsideGate(rxn1.InvMassExcEnergy,7.0,9.0)) // hump in exc. energy
    //     // { 
    //     ProcessFill();
    //     // } 
    // }

    //////////////////////////////////////////////
    // Final event reconstruction done here!
    // if(false) 
    if(true) 
    {
        // Select only monitoring ion chamber events
        // if(!InsideGate(IonCds.TRaw(),2490,2510.0) || !InsideGate(S1Trig.TRaw(),-1.0,1.0)) return;
        // if(!InsideGate(ic.ESeg.T(),119.0,123.0)) return;

        // Check anything inside silicon trigger OR anything outside IC down-scaled trigger
        if(InsideGate(IonCds.TRaw(),2490,2510)) return;
        if(!InsideGate(S1Trig.TRaw(),2100,2220)) return;

        // if(!InsideGate(S1Trig.T()-rf.TMod(), 547., 582.)) return;
        if(InsideGate(S1.Back.T(), 409., 413.) || InsideGate(S2.Back.T(), 409., 413.)) {} else { return; }

        if(!InsideGate(ic.ESeg.T(), 285., 360.)) return; 
        if(!InsideGate(ic_band1, ic.ESeg.ERaw(), ic.DESeg.ERaw())) return; 
        if(!InsideGate(protons, S2.Back.E()+S1.Back.E(), S2.Back.E())) return;

        ic.ReconstructCluster(ic.xgrid,ic.xcluster);
        ic.ReconstructCluster(ic.ygrid,ic.ycluster);
        if(ic.xcluster.Mult()>=2 || ic.ycluster.Mult()>=2) return;
        if(!InsideGate(ic.ycluster.E()-ic.xcluster.E(),591,1450)) return; 
        if(!InsideGate(ic.xcluster.E(),886,1482)) return; 
        
        ic.ReconstructHitPos(ic.xcluster.Ch(),ic.ycluster.Ch());

        // ic.xgrid.Print();
        // ic.xcluster.Print(); 
        // S2.FrontCluster.Print(); 

        if(S1.Back.Mult() < 1 || S1.Back.Mult() > 4) return;
        if(S1.Front.Mult() < 1 || S1.Front.Mult() > 4) return;
        if(S2.Back.Mult() < 1 || S2.Back.Mult() > 4) return;
        if(S2.Front.Mult() < 1 || S2.Front.Mult() > 4) return;

        TVector3 vS1 = S1.Ch2Vect(S1.Front.Ch(),S1.Back.Ch());
        TVector3 vS2 = S2.Ch2Vect(S2.Front.Ch(),S2.Back.Ch());
        if(!InsideGate(TMath::RadToDeg()*TMath::Abs(vS1.Theta()-vS2.Theta()),0,4.0)) return;  
        if(!InsideGate(TMath::RadToDeg()*TMath::Abs(vS1.Phi()-vS2.Phi()),0,45.0)) return; 

        S2.ReconstructClusters(S2); S1.ReconstructClusters(S1); 
        if(!InsideGate(TMath::RadToDeg()*TMath::Abs(vS1.Theta()-S1.Theta()),0,2.0)) return; 
        if(!InsideGate(TMath::RadToDeg()*TMath::Abs(vS1.Phi()-S1.Phi()),0,22.5)) return; 
        if(!InsideGate(TMath::RadToDeg()*TMath::Abs(S1.Theta()-S2.Theta()),0,2.0)) return;
        if(!InsideGate(TMath::Abs(S1.FrontCluster.E()-S1.BackCluster.E())/S1.BackCluster.E(),0,0.05)) return; 

        // if(!InsideGate(TMath::Abs(S2.FrontCluster.E()-S2.BackCluster.E())/S2.BackCluster.E(),0,0.05)) return; 
        // if(!InsideGate(TMath::Abs(S2.FrontCluster.E()-S2.BackCluster.E()),0,0.4)) return; 
        // if(!InsideGate(si_theta_corr, S2.Theta()*TMath::RadToDeg(), S1.Theta()*TMath::RadToDeg())) return; 
        // if(!InsideGate(si_phi_corr, S2.Phi()*TMath::RadToDeg(), S1.Phi()*TMath::RadToDeg())) return;    


        rxn1.DoPhysics(S1,S2,ic);
        inv_mass_ex = rxn1.InvMassExcEnergy;
        decay_sum = rxn1.Decay_Heavy.KE+rxn1.Decay_Light.KE;
        ProcessTree->Fill();
        inv_mass_ex = 0.0;

        ProcessFill();    
    }

    // S2.Print();
}// end of ProcessEvent





void Analyzer::ProcessFill(){ 
    if(saveTree) DataTree->Fill(); // Fill the TTree
    savedEvents+=1;

    /*Filling all Summary spectra like SpecTcl*/
    if(plot_summary_flag)
    {
        MyFill(table0,"Summary_S1_Front_ChRaw_vERaw",16,0,16,S1.Front.ChRaw(),4096,0,4096,S1.Front.ERaw());
        MyFill(table0,"Summary_S1_Back_ChRaw_vERaw",16,0,16,S1.Back.ChRaw(),4096,0,4096,S1.Back.ERaw());
        MyFill(table0,"Summary_S2_Front_ChRaw_vERaw",16,0,16,S2.Front.ChRaw(),4096,0,4096,S2.Front.ERaw());
        MyFill(table0,"Summary_S2_Back_ChRaw_vERaw",16,0,16,S2.Back.ChRaw(),4096,0,4096,S2.Back.ERaw());
        MyFill(table0,"Summary_S1_Back_ChRaw_vTRaw",16,0,16,S1.Back.ChRaw(),4096,0,4096,S1.Back.TRaw());
        MyFill(table0,"Summary_S1_Back_ChRaw_vT",16,0,16,S1.Back.ChRaw(),2000,0,1000,S1.Back.T());
        MyFill(table0,"Summary_S2_Back_ChRaw_vTRaw",16,0,16,S2.Back.ChRaw(),4096,0,4096,S2.Back.TRaw());
        MyFill(table0,"Summary_S2_Back_ChRaw_vT",16,0,16,S2.Back.ChRaw(),2000,0,1000,S2.Back.T());

        MyFill(table0,"Summary_S2_Back_E_vT",400,0,20,S2.Back.E(),2000,0,1000,S2.Back.T());

        MyFill(table0,"Summary_S2_Back_ChRaw_vTRaw",16,0,16,S2.Back.ChRaw(),4096,0,4096,S2.Back.TRaw());
        MyFill(table0,"Summary_ic_xgrid_ChRaw_vERaw",32,0,32,ic.xgrid.ChRaw(),4096,0,4096,ic.xgrid.E());
        MyFill(table0,"Summary_ic_ygrid_ChRaw_vERaw",32,0,32,ic.ygrid.ChRaw(),4096,0,4096,ic.ygrid.E());
    }

    if(plot_run_flag)
    {
        MyFill(table0,"run_ic_ESeg_ERaw",200,3200,3400,runNum,4096,0,4096,ic.ESeg.ERaw());
        MyFill(table0,"run_ic_DESeg_ERaw",200,3200,3400,runNum,4096,0,4096,ic.DESeg.ERaw());
        MyFill(table0,"run_ic_ESeg_TRaw",200,3200,3400,runNum,4096,0,4096,ic.ESeg.TRaw());
        MyFill(table0,"run_ic_E",200,3200,3400,runNum,200,0,100,ic.ESeg.E());
        MyFill(table0,"run_ic_DE",200,3200,3400,runNum,200,0,100,ic.DESeg.E());
        MyFill(table0,"run_rf",200,3200,3400,runNum,5000,0,500,rf.T());
    }
    /*Filling all Trigger/Timing related histograms*/
    if(plot_timing_flag)
    {
        MyFill(table0,"tdc1_6",2048,-1024,4096,tdc1[6]);
        MyFill(table0,"t_S1Trig-rfMod",2000,0,1000,S1Trig.T()-rf.TMod());
        MyFill(table0,"t_S1Trig-rf",2000,0,1000,S1Trig.T()-rf.T());
        MyFill(table0,"t_S1Trig-rf_vEt",2000,0,1000,S1Trig.T()-rf.TMod(),1000,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_rf_TRel_S2_vEt",2000,0,100,rf.TRel(S2.Back.T()),1000,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_rf",5000,0,500,rf.T());
        MyFill(table0,"t_rf_mod",5000,0,500,rf.TMod());
        MyFill(table0,"t_S1Trig",2054,-1024,4096,S1Trig.TRaw());
        MyFill(table0,"t_MCP",2054,-1024,4096,MCP.TRaw());
        MyFill(table0,"t_IonCds",2054,-1024,4096,IonCds.TRaw());
        MyFill(table0,"t_S1-rfMod",2000,0,1000,S1.Back.T()-rf.TMod());
        MyFill(table0,"t_S2-rfMod",2000,0,1000,S2.Back.T()-rf.TMod());
        MyFill(table0,"t_S1-rf",2000,0,1000,S1.Back.T()-rf.T());
        MyFill(table0,"t_S2-rf",2000,0,1000,S2.Back.T()-rf.T());
        MyFill(table0,"t_S1",2000,0,1000,S1.Back.T());
        MyFill(table0,"t_S2",2000,0,1000,S2.Back.T());
        MyFill(table0,"t_rf_TRel_S1",2000,0,100,rf.TRel(S1.Back.T()));
        MyFill(table0,"t_rf_TRel_S2",2000,0,100,rf.TRel(S2.Back.T()));
        MyFill(table0,"t_S1_Raw",2048,0,4096,S1.Back.TRaw());
        MyFill(table0,"t_S2_Raw",2048,0,4096,S2.Back.TRaw());
        MyFill(table0,"t_ic_T_vEt",1000,0,1000,ic.ESeg.T(),1000,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_IonC_Raw",4096,-1024,4096,ic.ESeg.TRaw());
        MyFill(table0,"t_IonC",1000,0,1000,ic.ESeg.T());
        MyFill(table0,"t_rf_rel_ic",2000,0,100,rf.TRel(ic.ESeg.T()));
        MyFill(table0,"t_rf_vEt",1000,0,500,rf.T(),500,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_rf_rel_ic_vEt",200,0,100,rf.TRel(ic.ESeg.T()),500,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_rf_mod_ic_vEt",5000,0,500,rf.TMod(),500,0,100,ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"t_IonC_vrf",500,0,1000,ic.ESeg.T(),5000,0,500,rf.T());
        MyFill(table0,"t_rf_MCP",500,-500,500,MCP.TRaw()-rf.TRaw());
        MyFill(table0,"t_rf_rel_MCP",500,-500,500,rf.TRel(MCP.T()));
    }

    /*Filling all SILICON related histograms */
    if(plot_si_flag)
    {
        MyFill(table0,"S1_Front_ERaw",2048,0,4096,S1.Front.ERaw());
        MyFill(table0,"S1_Front_Mult",16,0,16,S1.Front.Mult());
        MyFill(table0,"S1_Back_ERaw",2048,0,4096,S1.Back.ERaw());
        MyFill(table0,"S1_Back_Mult",16,0,16,S1.Back.Mult());

        MyFill(table0,"S2_Front_ERaw",2048,0,4096,S2.Front.ERaw());
        MyFill(table0,"S2_Front_Mult",16,0,16,S2.Front.Mult());
        MyFill(table0,"S2_Back_ERaw",2048,0,4096,S2.Back.ERaw());
        MyFill(table0,"S2_Back_Mult",16,0,16,S2.Back.Mult());

        MyFill(table0,"S1_Front_ChRaw",16,0,16,S1.Front.ChRaw());
        MyFill(table0,"S1_Back_ChRaw",16,0,16,S1.Back.ChRaw());
        MyFill(table0,"S2_Front_ChRaw",16,0,16,S2.Front.ChRaw());
        MyFill(table0,"S2_Back_ChRaw",16,0,16,S2.Back.ChRaw());
        MyFill(table0,"S1_Front_E",400,0,40,S1.Front.E());
        MyFill(table0,"S1_Back_E",400,0,40,S1.Back.E());
        MyFill(table0,"S2_Front_E",200,0,20,S2.Front.E());
        MyFill(table0,"S2_Back_E",200,0,20,S2.Back.E());
        MyFill(table0,"S2_FrontCluster_E",200,0,20,S2.FrontCluster.E());
        MyFill(table0,"S1_FrontCluster_E",400,0,40,S1.FrontCluster.E());
        MyFill(table0,"S1_ERaw_corr",1024,0,4096,S1.Front.ERaw(), 1024,0,4096,S1.Back.ERaw());
        MyFill(table0,"S2_ERaw_corr",1024,0,4096,S2.Front.ERaw(), 1024,0,4096,S2.Back.ERaw());
        MyFill(table0,"S1_EGain_corr",500,0,20,S1.Front.E_LocalGain(),500,0,20,S1.Back.E_LocalGain());
        MyFill(table0,"S2_EGain_corr",250,0,10,S2.Front.E_LocalGain(),250,0,10,S2.Back.E_LocalGain());
        MyFill(table0,"S1_E_corr",400,0,40,S1.Front.E(), 400,0,40,S1.Back.E());
        MyFill(table0,"S2_E_corr",400,0,20,S2.Front.E(), 400,0,20,S2.Back.E());
        MyFill(table0,"S1_FBdiff_BackE",400,0,40,S1.Back.E(),400,-5,5,S1.Front.E()-S1.Back.E());
        MyFill(table0,"S2_FBdiff_BackE",400,0,20,S2.Back.E(),400,-5,5,S2.Front.E()-S2.Back.E());
        MyFill(table0,"S1_FBdiff_FrontE",400,0,40,S1.Front.E(),400,-5,5,S1.Front.E()-S1.Back.E());
        MyFill(table0,"S2_FBdiff_FrontE",400,0,20,S2.Front.E(),400,-5,5,S2.Front.E()-S2.Back.E());

        MyFill(table0,"S1_Cluster_FBdiff_BackE",400,0,40,S1.BackCluster.E(),400,-5,5,S1.FrontCluster.E()-S1.BackCluster.E());
        MyFill(table0,"S2_Cluster_FBdiff_BackE",400,0,20,S2.BackCluster.E(),400,-5,5,S2.FrontCluster.E()-S2.BackCluster.E());
        MyFill(table0,"S1_Cluster_FBdiff_FrontE",400,0,40,S1.FrontCluster.E(),400,-5,5,S1.FrontCluster.E()-S1.BackCluster.E());
        MyFill(table0,"S2_Cluster_FBdiff_FrontE",400,0,20,S2.FrontCluster.E(),400,-5,5,S2.FrontCluster.E()-S2.BackCluster.E());

        MyFill(table0,"S1_Cluster_corr",500,0,20,S1.FrontCluster.E(), 500,0,20,S1.BackCluster.E());
        MyFill(table0,"S2_Cluster_corr",500,0,10,S2.FrontCluster.E(), 500,0,10,S2.BackCluster.E());
        MyFill(table0,"Si_Front_corr",16,0,16,S2.Front.Ch(),16,0,16,S1.Front.Ch());
        MyFill(table0,"Si_Front_corr_raw",16,0,16,S2.Front.ChRaw(),16,0,16,S1.Front.ChRaw());
        MyFill(table0,"Si_Back_corr",16,0,16,S2.Back.Ch(),16,0,16,S1.Back.Ch());

        MyFill(table0,"Si_PID_Raw",2048,0,8192,S2.Back.ERaw()+S1.Back.ERaw(), 1024,0,4096,S2.Back.ERaw());
        MyFill(table0,"Si_PID_Back",1200,0,60,S2.Back.E()+S1.Back.E(), 600,0,20,S2.Back.E());
        MyFill(table0,"Si_PID_Front",1200,0,60,S2.Front.E()+S1.Front.E(), 600,0,20,S2.Front.E());
        MyFill(table0,"Si_PID_Match",1200,0,60,S2.E()+S1.E(), 600,0,20,S2.E());
        MyFill(table0,"Si_PID_Cluster",1200,0,60,S2.BackCluster.E()+S1.BackCluster.E(), 600,0,20,S2.BackCluster.E());
        MyFill(table0,"Si_EvS1Theta",180,0,45,S1.Theta()*TMath::RadToDeg(), 400,0,40,S1.Back.E()+S2.Back.E());
        MyFill(table0,"Si_EvS2Theta",180,0,45,S2.Theta()*TMath::RadToDeg(), 400,0,40,S1.Back.E()+S2.Back.E());
        MyFill(table0,"S1_E_vRFmod",1000,0,200,rf.TMod(), 60,0,20,S1.Back.E());
        MyFill(table0,"S2_E_vRFmod",1000,0,200,rf.TMod(), 60,0,10,S2.Back.E());
        MyFill(table0,"S2_Theta", 180,0,45, S2.Theta()*TMath::RadToDeg());
        MyFill(table0,"S1_Theta", 180,0,45, S1.Theta()*TMath::RadToDeg());
        MyFill(table0,"S2_Front_Ch_vS1_Back_Ch",16,0,16,S2.Front.Ch(),16,0,16,S1.Back.Ch());
        MyFill(table0,"S1_Front_Ch_vS2_Back_Ch",16,0,16,S1.Front.Ch(),16,0,16,S2.Back.Ch());

        MyFill(table0,"si_theta_corr",180,0,45,S2.Ch2Vect(S2.Front.Ch(),S2.Back.Ch()).Theta()*TMath::RadToDeg(), 180,0,45,S1.Ch2Vect(S1.Front.Ch(),S1.Back.Ch()).Theta()*TMath::RadToDeg());
        MyFill(table0,"si_phi_corr",400,-200,200,S2.Ch2Vect(S2.Front.Ch(),S2.Back.Ch()).Phi()*TMath::RadToDeg(), 400,-200,200,S1.Ch2Vect(S1.Front.Ch(),S1.Back.Ch()).Phi()*TMath::RadToDeg());

        MyFill(table0,"si_rec_theta_corr",180,0,45,S2.Theta()*TMath::RadToDeg(), 180,0,45,S1.Theta()*TMath::RadToDeg());
        MyFill(table0,"si_rec_phi_corr",400,-200,200,S2.Phi()*TMath::RadToDeg(), 400,-200,200,S1.Phi()*TMath::RadToDeg());

        MyFill(table0,"s1_theta_corr",180,0,45,S1.Ch2Vect(S1.Front.Ch(),S1.Back.Ch()).Theta()*TMath::RadToDeg(), 
                180,0,45,S1.Theta()*TMath::RadToDeg());

        MyFill(table0,"s2_theta_corr",180,0,45,S2.Ch2Vect(S2.Front.Ch(),S2.Back.Ch()).Theta()*TMath::RadToDeg(), 
                180,0,45,S2.Theta()*TMath::RadToDeg());

        MyFill(table0,"s1_phi_corr",400,-200,200,S1.Ch2Vect(S1.Front.Ch(),S1.Back.Ch()).Phi()*TMath::RadToDeg(), 
                400,-200,200,S1.Phi()*TMath::RadToDeg());

        MyFill(table0,"s2_phi_corr",400,-200,200,S2.Ch2Vect(S2.Front.Ch(),S2.Back.Ch()).Phi()*TMath::RadToDeg(), 
                400,-200,200,S2.Phi()*TMath::RadToDeg());

        MyFill(table0,"S2_Pos",300,-75,75,S2.X(),300,-75,75,S2.Y());
        MyFill(table0,"S1_Pos",300,-75,75,S1.X(),300,-75,75,S1.Y());
        MyFill(table0,"S1_Radius",400,0,40,TMath::Sqrt(S1.X()*S1.X() + S1.Y()*S1.Y()));
        MyFill(table0,"S2_Radius",400,0,40,TMath::Sqrt(S2.X()*S2.X() + S2.Y()*S2.Y()));

        // for(int i=0;i<S1.Front.Mult();i++){
        //   MyFill(table0,Form("S1_Front_%d",(int)S1.Front.Ch(0)),16,0,16,S1.Front.Ch(i),4000,0,20,S1.Front.E(i));
        // }
        // for(int i=0;i<S2.Front.Mult();i++){
        //   MyFill(table0,Form("S2_Front_%d",(int)S2.Front.Ch(0)),16,0,16,S2.Front.Ch(i),2000,0,10,S2.Front.E(i));
        // }
        // for(int i=0;i<S1.Back.Mult();i++){
        //   MyFill(table0,Form("S1_Back_%d",(int)S1.Back.Ch(0)),16,0,16,S1.Back.Ch(i),4000,0,20,S1.Back.E(i));
        // }
        // for(int i=0;i<S2.Back.Mult();i++){
        //   MyFill(table0,Form("S2_Back_%d",(int)S2.Back.Ch(0)),16,0,16,S2.Back.Ch(i),2000,0,10,S2.Back.E(i));
        // }
    } // end of silicon plotting


    /* Filling NEUTRON related histrograms */
    if(plot_neutron_flag)
    {
        MyFill(table0,"Neut_ERaw",1024,0,4096,neut.det.ERaw());
        MyFill(table0,"Neut_dERaw",1024,0,4096,neut.det.dERaw());
        MyFill(table0,"Neut_PSD",300,0,3,neut.det.dERaw()/neut.det.ERaw(),1024,0,4096,neut.det.ERaw());
        MyFill(table0,"Neut_TRaw",1024,0,4096,neut.det.TRaw());
        MyFill(table0,"Neut_ChRaw",12,0,12,neut.det.ChRaw());
        MyFill(table0,"Neut_TRaw_Summary",12,0,12,neut.det.ChRaw(),1024,0,4096,neut.det.TRaw());
    }


    /* Filling Ion Chamber related histrograms */
    if(plot_ic_flag)
    {
        // for(int i=0;i<ic.xgrid.Mult();i++){
        //   MyFill(table0,Form("ic_xgrid_%d",(int)ic.xgrid.Ch(0)),16,0,16,ic.xgrid.Ch(i),512,0,4096,ic.xgrid.E(i));
        //   MyFill(table0,Form("ic_xgrid_%d_v%d",(int)ic.xgrid.Ch(i),(int)ic.xgrid.Ch(i+1)),512,0,4096,ic.xgrid.E(i),512,0,4096,ic.xgrid.E(i+1));
        // }
        // MyFill(table0,Form("ic_xgrid%d_ERaw_p1",(int)ic.xgrid.Ch(0)),32,0,32,ic.xgrid.Ch(0),512,0,4096,ic.xgrid.ERaw());

        MyFill(table0,"ic_xgrid_ERaw_summ",32,0,32,ic.xgrid.Ch(),512,0,4096,ic.xgrid.E_LocalGain());
        MyFill(table0,"ic_ygrid_ERaw_summ",32,0,32,ic.xgrid.Ch(),512,0,4096,ic.ygrid.E_LocalGain());

        MyFill(table0,"ic_xClust_yClust",512,0,4096,ic.xcluster.E(),512,0,4096,ic.ycluster.E());
        MyFill(table0,"ic_xClust1_xClust2",512,0,4096,ic.xcluster.E(0),512,0,4096,ic.xcluster.E(1));
        // MyFill(table0,Form("ic_x%d_y",(int)ic.xgrid.Ch()),512,0,4096,ic.xgrid.E_LocalGain(),512,0,4096,ic.ygrid.E_LocalGain());
        // MyFill(table0,Form("ic_x_y%d",(int)ic.ygrid.Ch()),512,0,4096,ic.xgrid.E_LocalGain(),512,0,4096,ic.ygrid.E_LocalGain());

        if(ic.xgrid.Mult()==1 && ic.ygrid.Mult()==1){
            MyFill(table0,"ic_GridGainMatch_m0",512,0,4096,ic.xgrid.E(),512,0,4096,ic.ygrid.E());
        }
        MyFill(table0,"ic_ESeg_ERaw",4096,0,4096,ic.ESeg.ERaw());
        MyFill(table0,"ic_DEDeg_ERaw",4096,0,4096,ic.DESeg.ERaw());
        MyFill(table0,"ic_ESeg_E",1000,0,100,ic.ESeg.E());
        MyFill(table0,"ic_DESeg_E",500,0,50,ic.DESeg.E());

        MyFill(table0,"ic_xGrid_mult",32,0,32,ic.xgrid.Mult());
        MyFill(table0,"ic_yGrid_mult",32,0,32,ic.ygrid.Mult());
        MyFill(table0,"ic_xClust_mult",32,0,32,ic.xcluster.Mult());
        MyFill(table0,"ic_yClust_mult",32,0,32,ic.ycluster.Mult());

        MyFill(table0,"ic_xgrid_ERaw",4096,0,4096,ic.xgrid.ERaw());
        MyFill(table0,"ic_ygrid_ERaw",4096,0,4096,ic.ygrid.ERaw());
        MyFill(table0,"ic_xGridE_yGridE",512,0,4096,ic.xgrid.E(),512,0,4096,ic.ygrid.E());
        MyFill(table0,"ic_xGridE_dESeg",512,0,4096,ic.xgrid.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"ic_xGridE_ESeg",512,0,4096,ic.xgrid.E(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"ic_yGridE_dESeg",512,0,4096,ic.ygrid.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"ic_yGridE_ESeg",512,0,4096,ic.ygrid.E(),1000,0,100,ic.ESeg.E());

        MyFill(table0,"ic_xClustE",1024,0,4096,ic.xcluster.E());  
        MyFill(table0,"ic_yClustE",1024,0,4096,ic.ycluster.E());
        MyFill(table0,"ic_xClustE_yClustE",512,0,4096,ic.xcluster.E(),512,0,4096,ic.ycluster.E());
        MyFill(table0,"ic_xClustE_dESeg",512,0,4096,ic.xcluster.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"ic_xClustE_ESeg",512,0,4096,ic.xcluster.E(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"ic_yClustE_dESeg",512,0,4096,ic.ycluster.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"ic_yClustE_ESeg",512,0,4096,ic.ycluster.E(),1000,0,100,ic.ESeg.E());

        MyFill(table0,"ic_xClustE_theta",1024,0,4096,ic.xcluster.E(),100,0,5,ic.Theta()*TMath::RadToDeg());
        MyFill(table0,"ic_yClustE_theta",1024,0,4096,ic.ycluster.E(),100,0,5,ic.Theta()*TMath::RadToDeg());
        MyFill(table0,"ic_xClustE_phi",1024,0,4096,ic.xcluster.E(),100,-200,200,ic.Phi()*TMath::RadToDeg());
        MyFill(table0,"ic_yClustE_phi",1024,0,4096,ic.ycluster.E(),100,-200,200,ic.Phi()*TMath::RadToDeg());

        MyFill(table0,"ic_PIDRaw",1024,0,2048,ic.ESeg.ERaw(),1024,0,2048,ic.DESeg.ERaw());
        MyFill(table0,"ic_PID_0",500,0,100,ic.ESeg.E(),250,0,50,ic.DESeg.E());
        MyFill(table0,"ic_PID_1",1200,0,120,ic.ESeg.E()+ic.DESeg.E(),500,0,50,ic.DESeg.E());
        MyFill(table0,"ic_PIDReco",1000,0,120,rxn1.ICReco(ic.ESeg.E()+ic.DESeg.E()),250,0,50,ic.DESeg.E());

        MyFill(table0,"ic_ESeg_Theta",100,0,5,ic.Theta()*TMath::RadToDeg(),1000,0,100,ic.ESeg.E());
        MyFill(table0,"ic_Theta",100,0,5,ic.Theta()*TMath::RadToDeg());
        MyFill(table0,"ic_Phi",100,-200,200,ic.Phi()*TMath::RadToDeg());
        MyFill(table0,"ic_Phi_vS1_Phi",20,-200,200,S1.Phi()*TMath::RadToDeg(),20,-200,200,ic.HitPos.Phi()*TMath::RadToDeg());
        MyFill(table0,"ic_Phi_vS2_Phi",20,-200,200,S2.Phi()*TMath::RadToDeg(),20,-200,200,ic.HitPos.Phi()*TMath::RadToDeg());
        MyFill(table0,"ic_xy_ChRaw",32,0,32,ic.xgrid.ChRaw(), 32,0,32,ic.ygrid.ChRaw());
        MyFill(table0,"ic_xy_Ch",32,0,32,ic.xgrid.Ch(), 32,0,32,ic.ygrid.Ch());
        MyFill(table0,"ic_xy_Cluster",32,0,32,ic.xcluster.Ch(), 32,0,32,ic.ycluster.Ch());
        MyFill(table0,"ic_xy",1000,-50,50,ic.HitPos.X(), 1000,-50,50,ic.HitPos.Y());
        MyFill(table0,"ic_E_vSiE",1000,0,100,ic.ESeg.E(),2000,0,20,S1.Back.E()+S2.Back.E());
        MyFill(table0,"ic_E_vS2E",1000,0,100,ic.ESeg.E(),2000,0,20,S2.Back.E());
        MyFill(table0,"ic_E_vS1E",1000,0,100,ic.ESeg.E(),2000,0,20,S1.Back.E());

        MyFill(table0,"ic_Reco",1000,0,120,rxn1.ICReco(ic.ESeg.E()+ic.DESeg.E()));
        MyFill(table0,"ic_Reco_plus_proton",60,70,120,S1.Back.E()+S2.Back.E()+rxn1.ICReco(ic.ESeg.E()+ic.DESeg.E()));
    }// end of ion chamber plotting


    /* Filling all rxn1 related histograms */
    if(plot_physics_flag)
    {// 5000keV/119bins=42keV/bin
        double Sp=5.514;
        int BW=120;

        MyFill(table0,"DecayQValueEst",BW,-1.,4.,rxn1.DecayQValueEst);
        MyFill(table0,"DecayQValueEst_vIonC_ESum",BW,-1.,4.,rxn1.DecayQValueEst,120,0,120,ic.ESeg.E()+ic.DESeg.E());

        MyFill(table0,"InvMassExcEnergy",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy);
        MyFill(table0,"InvMassExcEnergy_vS1Theta",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,50,0,30,S1.Theta()*TMath::RadToDeg());
        MyFill(table0,"InvMassExcEnergy_vS1Phi",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,50,-200,200,S1.Phi()*TMath::RadToDeg());
        MyFill(table0,"InvMassExcEnergy_vS2Theta",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,50,0,30,S2.Theta()*TMath::RadToDeg());
        MyFill(table0,"InvMassExcEnergy_vS2Phi",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,50,-200,200,S2.Phi()*TMath::RadToDeg());
        MyFill(table0,"InvMassExcEnergy_vEt",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,120,0,120, ic.ESeg.E()+ic.DESeg.E());
        MyFill(table0,"InvMassExcEnergy_vSiE",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,80,0,20, S1.E()+S2.E());
        MyFill(table0,"InvMassExcEnergy_vFrag",BW,Sp-1.,Sp+4.,rxn1.InvMassExcEnergy,40,70,110,rxn1.Fragment.KE);

        // MyFill(table0,"ic_Dep_vFragment",1200,0,120,ic.ESeg.E()+ic.DESeg.E(),70,70,140,rxn1.Fragment.KE);

        MyFill(table0,"ExcEnergy",BW,Sp-1.,Sp-4.,rxn1.DecayQValue+5.514);
        MyFill(table0,"DecayQValue",BW,-1,4,rxn1.DecayQValue);

        if(InsideGate(rxn1.InvMassExcEnergy,5.65,6.15)){ // only 3p state
            MyFill(table0,"Fragment_KE_3p",21,79,100, rxn1.Fragment.KE);
            MyFill(table0,"DecaySum_KE_3p",29,78,107, rxn1.Decay_Heavy.KE+rxn1.Decay_Light.KE);

            MyFill(table0,"decay_heavy_KE_3p",50,60,110,rxn1.Decay_Heavy.KE);
            MyFill(table0,"decay_light_KE_3p",200,0,20, rxn1.Decay_Light.KE);

            MyFill(table0,"DecayQValue_3p",120,-1,4,rxn1.DecayQValue);
            MyFill(table0,"Si_Theta_corr_3p",250,0,45,S2.Theta()*TMath::RadToDeg(), 225,0,45,S1.Theta()*TMath::RadToDeg());
            MyFill(table0,"Si_Sum_3p",500,0,20,S2.E()+S1.E());
            MyFill(table0,"ic_EdE_3p",500,0,100,ic.ESeg.E(),250,0,50,ic.DESeg.E());

            MyFill(table0,"ic_Esum_3p",100,0,100,ic.ESeg.E()+ic.DESeg.E());
            MyFill(table0,"Neut_Theta_3p",36,0,180,rxn1.Ejectile.LV.Theta()*TMath::RadToDeg());
            MyFill(table0,"Neut_ThetaCM_3p",36,0,180,rxn1.ThetaCM*TMath::RadToDeg()); 

            MyFill(table0,"ic_Dep_vFragment_3p",1200,0,120,ic.ESeg.E()+ic.DESeg.E(),70,70,140,rxn1.Fragment.KE);
        }
        if(InsideGate(rxn1.InvMassExcEnergy,6.1,6.46)){ // first 2p state
            MyFill(table0,"Fragment_KE_2p_1",40,70,110,rxn1.Fragment.KE);
            MyFill(table0,"DecaySum_KE_2p_1",40,70,110,rxn1.Decay_Light.KE+rxn1.Decay_Heavy.KE);
            MyFill(table0,"DecayQValue_2p_1",120,-1,4,rxn1.DecayQValue);

            MyFill(table0,"ic_Dep_vFragment_2p_1",1200,0,120,ic.ESeg.E()+ic.DESeg.E(),70,70,140,rxn1.Fragment.KE);
        }
        if(InsideGate(rxn1.InvMassExcEnergy,6.46,6.86)){ //only 3m state
            MyFill(table0,"Fragment_KE_3m",40,70,110, rxn1.Fragment.KE);
            MyFill(table0,"DecayQValue_3m",120,-1,4,rxn1.DecayQValue);

            MyFill(table0,"ic_Dep_vFragment_3m",1200,0,120,ic.ESeg.E()+ic.DESeg.E(),70,70,140,rxn1.Fragment.KE);
        }   
        if(InsideGate(rxn1.InvMassExcEnergy,7.2,7.66)){ // second 2p
            MyFill(table0,"Fragment_KE_2p_2",40,70,110, rxn1.Fragment.KE);
            MyFill(table0,"DecayQValue_2p_2",120,-1,4,rxn1.DecayQValue);
        }
    }



} // end of ProcessFill function



/******************************************************
  This is the big loop that goes over all the data!
  then grabs the data and fills the vector
 ******************************************************/
void Analyzer::Loop(){
    const ULong64_t nevents = data_chain->GetEntries();
    std::cout<<"Number of Events in root file are: "<<nevents<<std::endl;

    for(ULong64_t event=0; event < nevents; event++){
        data_chain->GetEntry(event);
        GetEventEntry();
        ProcessEvent();

        if(event%(nevents/100) == 0){
            std::cout<<100*event/nevents<<"% Completed -> "<<100.0*savedEvents/nevents<<"% Data Saved"<<std::flush<<"\r";
        }
    }
}


void Analyzer::Save(){ cout<<"\nAnalyzer Writing RootTree\n";
    out_file->cd();
    // DataTree->Write();
    table0->Write("",TObject::kOverwrite);
    std::cout<<"Number of events passed through gates: "<<savedEvents<<std::endl;
}




#endif
