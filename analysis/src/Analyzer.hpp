#ifndef __Analyzer__hpp
#define __Analyzer__hpp

#include "Detector.hpp"
#include "Physics.hpp"
#include "VariableMap.hpp"
#include "SiliconDetector.hpp"
#include "IonChamber.hpp"
#include "NeutronDetector.hpp"
#include "Timing.hpp"
#include "Global.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCutG.h"
#include "THashTable.h"

#include "catima/gwm_integrators.h"

using namespace std;

class Analyzer{
private:
  ULong64_t savedEvents = 0;
public:
  Analyzer(){};
  ~Analyzer(){};
  Analyzer(TChain*,TFile*,THashTable*);

  // Define root objects
  TTree *DataTree = new TTree("DataTree","DataTree");  
  // TTree *CleanTree = new TTree("CleanTree","CleanTree"); 
  TTree *ProcessTree = new TTree("ProcessTree","ProcessTree");  

  THashTable *table0 = new THashTable();
  TFile *out_file = new TFile();
  TChain *data_chain = new TChain();  

  // RxnTarget target, ic_window;
  // catima::Material target = catima::get_material(catima::material::CH2);
  catima::Material target = catima::Material({ {2,1,2},{12,6,1} });
  catima::Material kapton = catima::get_material(catima::material::Kapton);
  catima::Layers ion_chamber;
  catima::Material ic_dl = catima::get_material(catima::material::Isobutane);
  catima::Material ic_x = catima::get_material(catima::material::Isobutane);
  catima::Material ic_y = catima::get_material(catima::material::Isobutane);
  catima::Projectile beam1 = catima::Projectile(25,13,13);
  catima::Projectile beam2 = catima::Projectile(24,12,12);
  catima::Projectile beam3 = catima::Projectile(24,12,11);

  // Define analysis objected
  VariableMap VarMap;
  Physics rxn1;
  Float_t inv_mass_ex;
  Float_t decay_sum;

  SiliconDetector S2,S1; 
  IonChamber ic;
  NeutronDetector neut;
  Timing rf, MCP, IonCds, S1Trig;

  // Define all the graphical cuts
  TCutG *protons, *alphas;
  TCutG *si_theta_corr, *si_phi_corr;
  TCutG *ic_band1;
  TCutG *ds_al25_band;
  TCutG *ic_band2;
  TCutG *ic_cluster_corr;
  TCutG *ic_t_Et;
  TCutG *ic_p1, *ic_p2,*ic_p3;
  TCutG *InvMassExcEnergy_vFrag;

  // Plotting and running flags
  bool plot_si_flag=true;
  bool plot_ic_flag=true;
  bool plot_timing_flag=true;
  bool plot_physics_flag=true;
  bool plot_neutron_flag=true;
  bool plot_run_flag=true;
  bool plot_summary_flag=true;
  bool saveTree=true;

  //------------------------------------------------------------------------
  //Change what modules(arrays) and scalers(double) are read from DataFile
  //------------------------------------------------------------------------
  // Define any local Analyzer variables
  Int_t runNum=0;
  // Float_t RunNum[32];
  Float_t adc1[32];
  Float_t adc2[32];
  Float_t adc3[32];
  Float_t adc4[32];
  Float_t adc5[32];
  Float_t adc6[32]; //added for JBaker's data
  Float_t tdc1[32];
  Float_t tdc2[32];
  Float_t tdc3[32];
  Float_t qdc3[32];

  void LoadInputData();
  void LoadGates();
  void ReadGates(std::string);
  void SetCalibrations();
  void SetReactionInfo();
  void SetBranches();
  void ProcessEvent();
  void ProcessFill();
  void Loop();
  void Save();
  void GetEventEntry();

};

#endif
