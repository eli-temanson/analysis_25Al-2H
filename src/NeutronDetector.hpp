#ifndef __NeutronDetector__hpp
#define __NeutronDetector__hpp

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

// const static unsigned int kNEUTNODE(17); //node starts from 17 -> 32
// const static unsigned int kHOUSING(16);
// const static unsigned int kFHOUSING(15);
// const static unsigned int kTHOUSING(14);
// const static unsigned int kBHOUSING(50); // Node starts from 50 -> 65
// const static unsigned int kSILICONE(70); // Node starts from 70 -> 85
// const static unsigned int kABSORBER(5);
// const static unsigned int kORING(4);
// const static double kMidThick(13);
// const static double kFrontThick(101.72); //mm

class NeutronDetector{
 protected:
  std::string _Name;
  int _Num;

  double _T0;
  double _T1;
  double _TLin; 
  double _TShift;
  double _ShortOffset;
  int _APos;
  TVector3 _HitPos;
  TVector3 _Pos;
  double _PhePol2;
  double _PheLin;
  double _PheShift;
  double _ZRot;

  //NeutDetector Specific Parameterss 
  double _QLong;
  double _QShort;
  double _PSD;
  double _T;
  double _TRel; //relative to other detectors

 public:
  NeutronDetector();
  ~NeutronDetector();
  void init(std::string,const int&);
  std::string Name() const {return _Name;}
  int NumOfNeut() const {return _Num;}
  void Reset();
  void SetCalibrations(VariableMap& VarMap);

  Detector det,PosGrid;  
  
  void InsertPSDHit(const double& _ql, const double& _qs, const double& t = 0.0);
 
  double PSD() const ;
  double QLong() const ;
  double QShort() const;
  double QShortOffset() const;
  double QkeVee() const;  
  double T() const;
  double TLocal() const;
  double TRaw() const;
  double TRel() const;
  
  // double PosX()const;
  // double PosY()const;
  // double SumPosE()const;

  // double nKE(double tof) const;
  // double nKE_R(double tof) const;
  // double GetRadius()const ;
  // double GetThickness()const;
  // double GetThreshold()const;
  
  // double ZEst()const;
  // TVector3 GetHitPos()const;
  // int HitCounter() const ;
  // int GetArrayPos() const;
  // TVector3 GetPosVect() const; 

  // double Q_value_est(double tof, double m1, double m2,
	// 	                double beam_e, double& hi_KE, double& Q_val);

  // bool inDet(const TVector3& v , const TVector3& beamspot = TVector3(0,0,0));
  // double H_hit(TLorentzVector& inLV);
  // int NeutIn(TLorentzVector nLV,double& t,double& e, const TVector3& beamspot = TVector3(0,0,0));
  // int NeutInAfterScattering(TLorentzVector &nLV,double& t,double& e);
  // double C_hit(TLorentzVector& inLV);
  // void Reset();
  // void Build();
  // void SetCalibrations(VariableMap& detvar);
  // void SetThreshold(double thresh){ _Threshold = thresh;}
  // double CalculateTRel(const double &tfirst);
  // int IsANeutron();
  // int IsAGamma();
  // int HitID();
  // int SimCounter(){return fCounter;}
  // int SimCounterCarbon(){return fCounterCarbon;}

  // TGeoVolume *CreateCrystal(TGeoVolume *parent,const int& node, const TVector3 & shift = TVector3(0,0,0));
  // TGeoVolume *CreateCrystal(const int& node, const TVector3 & shift = TVector3(0,0,0));
  // TGeoVolume *CreateSilicone(TGeoVolume *parent,const int& node, const TVector3 & shift = TVector3(0,0,0));
  // TGeoVolume *CreateSilicone(const int& node, const TVector3 & shift = TVector3(0,0,0));


  // double QLong() const{return _QLong;}
  // double QShort() const{return _QShort;}
  // double QShortOffset() const{return (_QShort + _ShortOffset);}
  // double TRaw() const {return _T;}
  // double TLocal() const {return _T>0 ? ((_T * _T1) + _T0) : 0;}
  // double TRel() const {return _TRel;}
  // double PSD() const {return _PSD;}
  // double GetRadius()const {return _Radius;}
  // double GetThickness()const {return _Thickness;}
  // double GetThreshold()const{return _Threshold;}
  // TVector3 GetHitPos()const{return _HitPos;}
  // int HitCounter() const {return _Counter;}
  // int GetArrayPos() const{return _APos;}
  // TVector3 GetPosVect() const{return _Pos;}
  // double ZEst() const {return _Pos.Z() - (_Thickness/2);}


};


// class NeutronDetectorArray:public BaseClass{
// private:
// public:
//   double _Tfirst;
//   int _Detfirst;
//   int _Tmult;
//   int _Mult;

//   //TGeoVolume *gRNVolume;  

//   std::vector<TVector3> _Pos;
//   std::vector<double> _Qlong;
//   std::vector<double> _PSD;
//   std::vector<double> _T;
//   std::vector<int> _Detlist;

//   NeutronDetectorArray(const TString & name = "");
//   int ReconstructHits(NeutCollection& in);
//   int InsertHit(const double& q_long, const double& q_short, const double& q_T,
//                 const TVector3& _Pos,const int& index);
  
//   inline double T(int i=0) const {return i< _Mult ? _T[i] :0;}
 
//   int DetFirst() const {return _Detfirst;};
//   int Mult() const {return _Mult;}
//   int TMult() const {return _Tmult;}
//   int Det(int i=0) const {return i<_Mult ? _Detlist[i] : -1;}

//   virtual void Reset();

// };



  
// class AluminiumHousing{
// protected:
//   TVector3 _HitPos; //where on the front of the housing did the particle hit first
//   TVector3 _Pos; // detector position in world
//   int _AluminiumCounter; //how many aluminum interactions before termination
//   int _AbsorberCounter; //how many aluminum interactions before termination
//   double _ELost;
//   double _Outerrad; //outer dimensions of entire housing
//   double _Innerrad; //beam pipe
//   double _MinThreshold; //Will find the lowest n det threshold
//   double _Thickness; //thickness of housing
//   double _Threshold; //neutron threshold
//   double _Rotz;
//   int _HitID;

//   AngularDistribution angdis[4];



// public:
//   AluminiumHousing(const std::string& name = "");

//   int AluminiumCounter(){return fAluminiumCounter;}
//   int AbsorberCounter(){return fAbsorberCounter;}
//   bool inDet(const TVector3& v, const TVector3& beamspot);
//   int NeutIn(TLorentzVector nLV,double &t,double& e, const TVector3& beamspot);
//   double GetNextInteractionPoint(double nKE /*Lab in MeV*/);
//   double Al_hit(TLorentzVector&inLV);
//   double GetAbsorberInteractionPoint(double nKE /*Lab in MeV*/);  
//   double Absorber_hit(TLorentzVector&inLV);
  
//   TVector3 GetPosVect() const;
//   TVector3 GetHitPos() const;

//   TGeoVolume *CreateVolume(TGeoVolume*parent,const int& node, const TVector3&shift  = TVector3(0,0,0));
//   TGeoVolume *CreateVolume(TGeoVolume*parent);
//   TGeoVolume *CreateFrontVolume(TGeoVolume*parent);
//   TGeoVolume *CreateFrontVolume(TGeoVolume*parent,const int& node, const TVector3&shift  = TVector3(0,0,0));
//   TGeoVolume *CreateTargetVolume(TGeoVolume*parent);
//   TGeoVolume *CreateTargetVolume(TGeoVolume*parent,const int& node, const TVector3&shift  = TVector3(0,0,0));
//   TGeoVolume *CreateBackingVolume(TGeoVolume*parent);
//   TGeoVolume *CreateBackingVolume(TGeoVolume*parent,const int& node, const TVector3&shift  = TVector3(0,0,0));
//   TGeoVolume *CreateORING(TGeoVolume*parent,const int&node =kORING,const TVector3&shift = TVector3(0,0,0));


//   void Reset();
//   void SetCalibrations(VariableMap& detvar);

// };
  

// inline TVector3 AluminiumHousing::GetPosVect() const{return fPos;}
// inline TVector3 AluminiumHousing::GetHitPos() const{return fHitPos;}
// int ReconstructTREL(NeutCollection& in,int&t_mult,double&t_first,int& det_first);
// int PositionMap(int slot,TVector3 & pos, double rot = 0);
// TGeoVolume *CreateRNVolume();
// //  TGeoVolume *CreateAbsorbers();


// extern TGeoVolume * gRNVolume; 



// namespace sim{
//   extern double stepcounter;
//   extern double steptime;
//   extern double z_pos;
//   extern bool fgIsDirect;
//   extern bool fgNeutHit;
// }

#endif
