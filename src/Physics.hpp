#ifndef __Physics__hpp
#define __Physics__hpp

#include <cmath>

#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "VariableMap.hpp"
#include "Global.hpp"

#include "SiliconDetector.hpp"
#include "IonChamber.hpp"
#include "NeutronDetector.hpp"
#include "Timing.hpp"

// #include "catima/catima.h"
#include "catima/gwm_integrators.h"

struct Particle {
    double M,InvMass,KE,E,Ecm;
    double X,Y,Z;
    double Theta,Phi;
    double Px,Py,Pz,P;
    TLorentzVector LV;
    TVector3 Pos;
};

class Physics{
protected:
  double mass[6];
  double _p2,_p1,_p0;
  double FragEst=0.0;
  double Beam_Sigma=0.0,Beam_TKE=0.0;
  TRandom3 *Rndm = new TRandom3();

public:
  Physics();
  ~Physics(){}

  Particle Beam,Target,Ejectile,Fragment,Decay_Light,Decay_Heavy;
  VariableMap MassTable;

  // RxnTarget target, ic_window;
  // catima::Material target_half = catima::get_material(catima::material::CH2);
  catima::Material target_half = catima::Material({ {2,1,2},{12,6,1} });
  catima::Material kapton = catima::get_material(catima::material::Kapton);
  catima::Material s2 = catima::get_material(catima::material::SiO2);

  catima::Layers ion_chamber;
  catima::Material ic_dl = catima::get_material(catima::material::Isobutane);
  catima::Material ic_x = catima::get_material(catima::material::Isobutane);
  catima::Material ic_y = catima::get_material(catima::material::Isobutane);
  // catima::Material ic_de = catima::get_material(catima::material::Isobutane);
  // catima::Material ic_e = catima::get_material(catima::material::Isobutane);

  catima::Projectile light = catima::Projectile(1,1);
  catima::Projectile heavy = catima::Projectile(25,13);

  void SetReaction(const std::string& beam,const std::string& target,const std::string& ejectile, 
		              const std::string& fragment,const std::string& product,const std::string& heavy);

  void DoPhysics(const SiliconDetector &S1,const SiliconDetector &S2,const IonChamber &ic);
 
  void CalcQValue();
  void CalcDecayQValueEst();
  void CalcDecayQValue();
  void CalcExcEnergy();
  void CalcInvMassExcEnergy();
  void Calc_CM_Angle();
  void CalcRecoilQValue(const double&,const double&);
  void Print();

  double SeparationEnergy();

  double QValue = 0.0;
  double DecayQValueEst = 0.0;
  double DecayQValue = 0.0;
  double ExcEnergy = 0.0;
  double InvMassExcEnergy = 0.0;
  double RecoilQValue = 0.0;
  double ThetaCM = 0.0;
  double Sp=0.0;
  double gamma,term1,term2,term3;
  double Ecm,Pcm;

  void SetBeamKE(const double&);
  void SetBeamSigma(const double&);
  void SetFragEst(const double&);
  void SetICReco(const double&,const double&,const double&);
  double ICReco(const double&);

};


#endif
