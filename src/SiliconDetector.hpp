#ifndef __SiliconDetector__hpp
#define __SiliconDetector__hpp

//C and C++ libraries.
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

//ROOT libraties
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "Detector.hpp"
#include "VariableMap.hpp"

class SiliconDetector{
protected:
  std::string name;

  int frontNum, backNum;
  TVector3 normVect;
  TVector3 posVect;
  double rot[3] = {0.0, 0.0, 0.0};

  int ChFrontRejection[16], ChBackRejection[16]; 
  int RejectionNum;

  //For the old cluster class
  // double energyMatchBack=0.;
  // double energyMatchFront=0.;
  int frontMatchStat=0;
  int backMatchStat=0;
  bool UseBackEnergy=false;
  double matchEpsilon=0.;
  double matchDelta=1.;
  bool f_cluster=false,b_cluster=false;

  //From online data sheet in [mm]
  // double S2INNERRAD = 23.06/2.; 
  // double S2OUTERRAD = 70.0/2.;
  // double S1INNERRAD = 48.0/2.;
  // double S1OUTERRAD = 96.0/2.;

  double innerRad=1.;
  double outerRad=1.;
  double ringPitch=1.;
  double deltaPhi=1.;

  //From the old cluster class
  std::vector<double> Ch_back,Ch_front;
  std::vector<double> EMatch,TMatch;

  TRandom3 *Rndm = new TRandom3();


public:
  SiliconDetector(){}; //Constructor
  ~SiliconDetector(){}; //Destructor

  void init(std::string, const int&, const int&); 
  void Reset();
  void Print();
  std::string Name() const {return name;}

  Detector Front,Back,FrontCluster,BackCluster;
  
  std::vector<TVector3> Pos;
  TVector3 Ch2Vect(double,double);
  TVector3 GetPosVect() const {return posVect;}  
  void CalcNormVect();
  bool InDet(const TVector3& incident, const TVector3& beam = TVector3(0,0,0));
  void SetCalibrations(VariableMap& VarMap);
    
  void ReconstructClusters(SiliconDetector& Det);
  void CalcHitPos(){ Ch2Vect(Front.Ch(),Back.Ch()); }
 
  double E(int i=0) const {return EMatch[i];}
  double T(int i=0) const {return TMatch[i];}
  double Phi(int i=0) const {return Pos[i].Phi();}
  double Theta(int i=0) const {return Pos[i].Theta();}
  double X(int i=0) const {return Pos[i].X();}
  double Y(int i=0) const {return Pos[i].Y();}
  double Z(int i=0) const {return Pos[i].Z();}
  double ChFront(int i = 0) const {return Ch_front[i];}
  double ChBack(int i = 0) const {return Ch_back[i];}
  double MatchDelta() const{ return matchDelta;}
  double MatchEpsilon() const{ return matchEpsilon;}

  double ThetaMin();
  double ThetaMax();

};

class SiliconArray{
private:
  int NumOfSi;
  std::string name;
  double thickness;
  //EnergyLoss _ELossTableSiA;
  std::string PID;

public:  
  //SiArray(const std::string& name = "SiArray", int num = 2, const std::string& PID="1H");
  //void ReconstructHits(S2ClusterCollection& si_c_);
  SiliconArray(){}; //Constructor
  ~SiliconArray(){}; //Destructor
  void init(std::string); 

  Detector Si;

  std::string Name() const {return name;}
  void Reset();
  double ERecoAB();
  double ERecoA();
  
  void SetCalibrations(VariableMap&);
};






#endif
