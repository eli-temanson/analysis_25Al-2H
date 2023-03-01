#ifndef __Detector__hpp
#define __Detector__hpp

#include <iostream>
#include <vector>
#include <array>
// #include <bits/stdc++.h> 

#include "TROOT.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TFile.h"

#include "VariableMap.hpp"
#include "Global.hpp"

struct Set{
  double ERaw, TRaw, dERaw;
  double ChRaw;
  double E, T;
};
struct Calibration{
  double A0,A1,E0,E1,T0,T1;
  double TOffset, EOffset, EScale;
  int ChMap;
};

class Detector{
protected:
  int _num = 0, _mult = 0;
  std::string _name;
   
  double _LowELimit = -4096.;
  double _HighELimit = 4096.;
  double _ELin = 1.0;
  double _EShift = 0.0;
  double e_calib_p2 = 0.0;
  double _TLin = 1.0;
  double _TShift = 0.0;

  //Specifiy if your detector needs Run Dependent Calibration Default setting is false/No
  int _RDC = 0;
  double ResetValue = -500.;
     
public:
  Detector(){};
  ~Detector(){};
  
  void Initialize(std::string, const int&);
  void Reset();

  std::array<Set,32> Data {};
  std::array<Calibration,32> Calib {};
  
  int Mult() const { return _mult; } 
  
  double ERaw(int i=0) const;
  double dERaw(int i=0) const;
  double TRaw(int i=0) const;
  double ChRaw(int i=0) const;
  double Ch(int i=0) const;
  double E_Offset(int i=0) const;
  double E_LocalGain(int i=0) const; 
  double E_LocalCalib(int i=0) const;  
  double dE_LocalGain(int i=0) const;
  double dE_LocalCalib(int i=0) const;  
  double E(int i=0) const;
  double dE(int i=0) const;
  double T_Offset(int i=0) const;
  double T_LocalCalib(int i=0) const;
  double T(int i=0) const;
  
  void SetELin(const double& elin){ 
    _ELin=elin; 
  }
  void SetEShift(const double& eshift){ 
    _EShift=eshift; 
  }
  void SetELimits(const double& elow, const double& ehigh){
    _LowELimit=elow; _HighELimit=ehigh;
  }
  double EMin(){ return _LowELimit; }
  double EMax(){ return _HighELimit;}
    
  int NumOfCh() const{return _num;}
  std::string GetName() const{return _name;}
  
  void InsertHit(const double&, const double&,const double&);
  void InsertHit(const double&, const double&, const double&, const double&);
  void SortByEnergy();
  void SortByTime();
  void SortNeut();
  void SortByCh();
  
  void SetCalibrations(VariableMap& VarMap);
  void Print();

};







#endif
