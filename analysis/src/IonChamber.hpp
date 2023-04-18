/******************************************************************
This class to be modified using LSU methods for the ion chamber
******************************************************************/

#ifndef __IonChamber__hpp
#define __IonChamber__hpp

//C and C++ libraries.
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <memory>
#include <list>
#include <algorithm>

//ROOT libraties
#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TObject.h>
#include <TLine.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TPolyLine3D.h>
#include <TVector3.h>

#include "VariableMap.hpp"
#include "Detector.hpp"
#include "Global.hpp"

class IonChamber{
protected:
  std::string _Name;
  double _WindowMax = 37;
  double _ZPos;
  double _YPos;
  double _XPos;
  double _WireSeparation;
  TVector3 _RotVect;
  TVector3 _NormVect;
  std::vector<TVector3> _ClusterPos;
  double _IrisOpening;
  TRandom3 *ZeroToOne = new TRandom3();
  double _p2=0.0, _p1=0.0, _p0=0.0;

public:
  IonChamber();
  ~IonChamber();
  std::string Name() const {return _Name;}
  
  Detector ESeg;
  Detector DESeg;
  Detector xgrid;
  Detector ygrid;
  Detector xcluster;
  Detector ycluster;

  TVector3 HitPos;

  void init(std::string);
  void Reset();
  void SetCalibrations(VariableMap& VarMap);
  void Print();
  
  double SumE_X() const;//sum all wires xgrid section
  double SumE_Y() const;//sum all wires ygrid section
  double PosX();
  double PosY();
  double PosXAdj();
  double PosYAdj();
  void ReconstructXClusters();
  void ReconstructYClusters();
  void ReconstructHitPos(double chx,double chy);
  void ReconstructCluster(Detector& grid,Detector& cluster);

  TVector3 GetPosVect() const{return TVector3(_XPos,_YPos,_ZPos) ;}
  double Theta() const {return HitPos.Theta(); } //returns radians
  double Phi() const {return HitPos.Phi(); } //returns radians

  bool InDet(const TVector3& incident, const TVector3& beam = TVector3(0,0,0));
  bool Vect_to_ch(const TVector3&, double&, double&, const TVector3& = TVector3(0,0,0));
  void CalcNormVect();

  double E() const {return ESeg.E()+DESeg.E()+xgrid.E()+ygrid.E();}
  double DE() const  {return DESeg.E()+xgrid.E()+ygrid.E();}
  double ERaw() const {return ESeg.ERaw()+DESeg.ERaw()+ygrid.ERaw()+ygrid.ERaw();}
  double DERaw() const {return DESeg.ERaw()+ygrid.ERaw()+ygrid.ERaw();}
  double SumE_Pos()const {return SumE_X()+SumE_Y();} //sum-both pos grid sections
  double T() const  {return ESeg.T();}
  double TRaw() const {return ESeg.TRaw();}

  double GetIncomingEnergy() const {return _p2*(ESeg.E()+DESeg.E())*(ESeg.E()+DESeg.E()) + _p1*(ESeg.E()+DESeg.E()) + _p0;}

  double X() const {return HitPos.X();}
  double Y() const {return HitPos.Y();}
  double Z() const {return HitPos.Z();}


};


#endif
