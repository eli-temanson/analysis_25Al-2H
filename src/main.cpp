
#ifndef main_cpp
#define main_cpp

#include "Analyzer.hpp"

#include <iostream>
#include <filesystem>
#include <string>

#include "TObject.h"
#include "TROOT.h"
#include "TApplication.h"
#include <TChain.h>
#include "THashTable.h"
#include "TFile.h"

int main(int argc,char *argv[])
{
    std::cout<<"-------------------------------------------------\n";
    std::cout<<"--- Analysis Software for ResoNeut ---\n";
    std::cout<<"-------------------------------------------------\n"; 

    TApplication app("app",&argc,argv);

    // std::cout << "Current Path: " << std::filesystem::current_path() << std::endl;
    // std::filesystem::current_path("25Al+d_2014/analysis");
    // std::cout << "Setting Path: " << std::filesystem::current_path() << std::endl;

    // auto data_path = std::filesystem::path("/mnt/data0/2014_06_25Al_dn_jbaker/ds_output.root").string().c_str();

    auto out_file = new TFile("C:\\Users\\elite\\OneDrive\\Research\\analysis_25Al-2H-final\\data\\histograms\\output_2023_03_14.root","RECREATE");
    // auto out_file = new TFile("/mnt/data0/2014_06_25Al_dn_jbaker/output_2023_01_09.root","RECREATE");
    auto data_chain = new TChain("DataTree");
    auto hist_table = new THashTable();  

    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run322[3456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run323[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run324[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run325[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run326[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run327[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run328[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run329[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run330[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run331[0123456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run332[013456789].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run333[134567].root");

    // Test
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run322[3].root");
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run323[0123456789].root");

    // Used for Down-scaled 
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/run323[012].root");
    
    // Used for final results 
    // data_chain->Add("/mnt/data0/2014_06_25Al_dn_jbaker/Filtered_noDS.root");
    data_chain->Add("C:\\Users\\elite\\OneDrive\\Research\\analysis_25Al-2H-final\\data\\filtered\\Filtered_noDS.root");

    Analyzer b(data_chain, out_file, hist_table);
    b.LoadInputData();
    b.LoadGates();
    b.SetReactionInfo();
    b.SetCalibrations();
    b.SetBranches(); //Initalize variables used exclusively in Analysis.cpp
    b.Loop(); //Loop through the data_chain events and execute Process/ProcessFill
    b.Save();


    std::cout<<"Output File: "<<out_file->GetName()<<std::endl;
    out_file->Write(out_file->GetName(), TObject::kOverwrite);
    out_file->Close(); 

    // data_chain->Reset();

    return 0;

}

#endif
