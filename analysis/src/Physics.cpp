/*************************************************************
* Physics Calculations 
*
* 0 + 1 -> 2 + 3* -> 2 + (4+5)
* Beam + Target -> Ejectile + Fragment* -> Ejectile + (Decay_Light + Decay_Heavy)
* All Everything should be in MeV via natural units
*************************************************************/
#ifndef __Physics__cpp
#define __Physics__cpp

#include "Physics.hpp"

using namespace std;

/****************************/
  bool DEBUG = false; 
/****************************/
Physics::Physics()
{
  s2.thickness(0.014848);

  target_half.density(0.7987); // g/cm3
  target_half.thickness(0.000258); // g/cm2

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

  // ic_de.density(0.00011);
  // ic_de.thickness(0.00088);
  // ion_chamber.add(ic_de);

  // ic_e.density(0.00011);
  // ic_e.thickness(0.0022);
  // ion_chamber.add(ic_e);
}

void Physics::SetBeamKE(const double& e){ Beam_TKE = e;}
void Physics::SetBeamSigma(const double& esig){ Beam_Sigma = esig;}
void Physics::SetFragEst(const double& e){ FragEst = e; }
void Physics::SetICReco(const double& p2,const double& p1,const double& p0){ _p2=p2; _p1=p1; _p0=p0; }
double Physics::ICReco(const double& ic_e){ return _p2*ic_e*ic_e + _p1*ic_e + _p0; }

void Physics::DoPhysics(const SiliconDetector &S1,const SiliconDetector &S2,const IonChamber &ic){
  
  // std::cout<<"Stopping power, " << catima::dedx(light(S1.Back.E()+S2.Back.E()), target_half) << ", MeV/(g/cm2)" << std::endl;
  // Assign proton properties
  // Decay_Light.KE = S1.Back.E() + S2.Back.E();
  // Decay_Light.KE = (S1.Back.E() + S2.Back.E())*0.9923 + 0.0717;
  // Decay_Light.KE = (S1.BackCluster.E()+S2.BackCluster.E())*0.9923 + 0.0717;
  
  Decay_Light.Theta = S1.Theta();
  Decay_Light.Phi = S1.Phi();

  // double efinal = S1.Back.E() + S2.Back.E();
  double efinal = S1.BackCluster.E() + catima::reverse_integrate_energyloss(light(S1.BackCluster.E()/light.A), s2);
  
  // set default thickness
  target_half.thickness(0.000258); // g/cm2
  kapton.thickness(0.000994);
  ic_dl.thickness(0.00036);
  ic_x.thickness(0.00044);
  ic_y.thickness(0.00044);

  target_half.thickness(0.000258/cos(Decay_Light.Theta)); // g/cm2
  efinal = efinal + catima::reverse_integrate_energyloss(light(efinal/light.A), target_half);
  
  // reset half-target thickness for the ion chamber
  target_half.thickness(0.000258); // g/cm2

  Decay_Light.KE = efinal;

  // Decay_Light.KE = S1.Back.E() + S2.Back.E();
  Decay_Light.E = Decay_Light.KE + Decay_Light.M;
  Decay_Light.P = sqrt(Decay_Light.E*Decay_Light.E - Decay_Light.M*Decay_Light.M);
  Decay_Light.Px = Decay_Light.P*cos(Decay_Light.Phi)*sin(Decay_Light.Theta);
  Decay_Light.Py = Decay_Light.P*sin(Decay_Light.Phi)*sin(Decay_Light.Theta);
  Decay_Light.Pz = Decay_Light.P*cos(Decay_Light.Theta);
  Decay_Light.LV.SetPxPyPzE(Decay_Light.Px, Decay_Light.Py, Decay_Light.Pz, Decay_Light.E);

  // Assign heavy ion chamber properties
  // Decay_Heavy.KE = ICReco(ic.ESeg.E()+ic.DESeg.E());
  Decay_Heavy.Theta = ic.HitPos.Theta();  
  Decay_Heavy.Phi = ic.HitPos.Phi();
  // Decay_Heavy.Theta = 0.0;
  // Decay_Heavy.Phi = 0.0;

  // target_half.thickness(0.000258/cos(Decay_Heavy.Theta)); // g/cm2
  // kapton.thickness(0.000994/cos(Decay_Heavy.Theta));
  // ic_dl.thickness(0.00036/cos(Decay_Heavy.Theta));
  // ic_x.thickness(0.00044/cos(Decay_Heavy.Theta));
  // ic_y.thickness(0.00044/cos(Decay_Heavy.Theta));

  double ein = 0.0;
  ein = ic.DESeg.E() + ic.ESeg.E();
  ein = ein + catima::reverse_integrate_energyloss(heavy(ein/heavy.A), ic_y);
  ein = ein + catima::reverse_integrate_energyloss(heavy(ein/heavy.A), ic_x);
  ein = ein + catima::reverse_integrate_energyloss(heavy(ein/heavy.A), ic_dl);
  ein = ein + catima::reverse_integrate_energyloss(heavy(ein/heavy.A), kapton);
  ein = ein + catima::reverse_integrate_energyloss(heavy(ein/heavy.A), target_half);

  Decay_Heavy.KE = ein;

  Decay_Heavy.E = Decay_Heavy.KE + Decay_Heavy.M;
  Decay_Heavy.P = sqrt(Decay_Heavy.E*Decay_Heavy.E - Decay_Heavy.M*Decay_Heavy.M);
  Decay_Heavy.Px = Decay_Heavy.P*cos(Decay_Heavy.Phi)*sin(Decay_Heavy.Theta);
  Decay_Heavy.Py = Decay_Heavy.P*sin(Decay_Heavy.Phi)*sin(Decay_Heavy.Theta);
  Decay_Heavy.Pz = Decay_Heavy.P*cos(Decay_Heavy.Theta);
  Decay_Heavy.LV.SetPxPyPzE(Decay_Heavy.Px, Decay_Heavy.Py, Decay_Heavy.Pz, Decay_Heavy.E);

  // Assume the Beam is on axis
  Beam.KE = Rndm->Gaus(Beam_TKE, Beam_Sigma);
  Beam.E = Beam.KE + Beam.M;
  Beam.P = sqrt(Beam.E*Beam.E - Beam.M*Beam.M);
  Beam.LV.SetPxPyPzE(0,0,Beam.P,Beam.E);
  //Assign Stationary target
  Target.LV.SetPxPyPzE(0,0,0,Target.M);
 
  //Method 1 (invarient mass method)
  Fragment.LV = Decay_Light.LV + Decay_Heavy.LV;
  Fragment.E = sqrt(Fragment.LV.P()*Fragment.LV.P() + Fragment.M*Fragment.M);
  Fragment.KE = Fragment.E - Fragment.M;

  //Method 2 (Calculate Neutron Kinematics)
  // Ejectile.Px = Beam.LV.Px() - Decay_Light.Px;
  // Ejectile.Py = Beam.LV.Py() - Decay_Light.Py;
  // Ejectile.Pz = Beam.LV.Pz() - Decay_Light.Pz;
  // Ejectile.P = sqrt(Ejectile.Px*Ejectile.Px + Ejectile.Py*Ejectile.Py + Ejectile.Pz*Ejectile.Pz);
  // Ejectile.E = sqrt(Ejectile.P*Ejectile.P + Ejectile.M*Ejectile.M);
  // Ejectile.LV.SetPxPyPzE(Ejectile.Px,Ejectile.Py,Ejectile.Pz,Ejectile.E);
  // Ejectile.KE = Ejectile.E - Ejectile.M;

  // //Method 3 ()
  // boostvector = (Beam.LV + Target.LV).BoostVector();
  // Ejectile.LV =  - Fragment.LV;
  // Ejectile.E = sqrt(Ejectile.LV.P()*Ejectile.LV.P() + Ejectile.M*Ejectile.M);
  // Ejectile.KE = Ejectile.E - Ejectile.M;  

  CalcDecayQValueEst();
  CalcDecayQValue();
  CalcQValue();
  CalcInvMassExcEnergy();
  CalcExcEnergy();
  Calc_CM_Angle();

  if(DEBUG)
  {
    // std::cout<<"\n\n------------- DEBUG: Physics Variables ----------------\n";

  }
}


//Set Primary Reaction and automatically add one decay channel
void Physics::SetReaction(const string &beam,const string &target,const string &ejectile, 
  			                  const string &fragment,const string &light,const string &heavy){
  
  cout<<"Setting Reaction \n";
  if(light == "p") { cout<<"Error!! use 1H not p \n"; exit(0); }
  MassTable.LoadParams("assets/input/MassTable.in");
  MassTable.GetParam(beam, Beam.M); 
  MassTable.GetParam(target, Target.M); 
  MassTable.GetParam(ejectile, Ejectile.M);
  MassTable.GetParam(fragment, Fragment.InvMass); 
  MassTable.GetParam(light, Decay_Light.M); 
  MassTable.GetParam(heavy, Decay_Heavy.M);
  Fragment.M = Fragment.InvMass;
}

double Physics::SeparationEnergy(){ 
  Sp = Decay_Heavy.M + Decay_Light.M - Fragment.M; 
  return Sp;
}

void Physics::Calc_CM_Angle(){
  TVector3 boostvector;

  boostvector = (Beam.LV + Target.LV).BoostVector();
  Beam.LV.Boost(-boostvector); // boost into cm frame
  Target.LV.Boost(-boostvector);
  Fragment.LV.Boost(-boostvector);

  Ejectile.LV = Beam.LV + Target.LV - Fragment.LV;
  // Ejectile.E = sqrt(Ejectile.LV.P()*Ejectile.LV.P() + Ejectile.M*Ejectile.M);
  // Ejectile.KE = Ejectile.E - Ejectile.M;  

  ThetaCM = Ejectile.LV.Theta();

  Beam.LV.Boost(boostvector); // boost back into lab frame
  Target.LV.Boost(boostvector);
  Fragment.LV.Boost(boostvector);
  Ejectile.LV.Boost(boostvector);
}

void Physics::CalcDecayQValueEst(){
  double a = Decay_Light.KE*(1 + Decay_Light.M/Decay_Heavy.M);
  double b = FragEst*(1 - Fragment.M/Decay_Heavy.M);
  double c = 2.0/Decay_Heavy.M;
  double d = sqrt(Decay_Light.KE * Decay_Light.M * FragEst * Fragment.M); 

  DecayQValueEst = a - b - c*d*cos( Decay_Light.Theta );
}

void Physics::CalcDecayQValue(){
  double a = Decay_Light.KE*(1 + Decay_Light.M/Decay_Heavy.M);
  double b = Fragment.KE*(1 - Fragment.M/Decay_Heavy.M);
  double c = 2.0/Decay_Heavy.M;
  double d = sqrt(Decay_Light.KE * Decay_Light.M * Fragment.KE * Fragment.M); 

  DecayQValue = a - b - c*d*cos( Decay_Light.Theta);
}

void Physics::CalcQValue(){
  QValue = Decay_Light.KE + Decay_Heavy.KE + Ejectile.KE - Beam.KE - Target.M;
}

void Physics::CalcExcEnergy(){
  ExcEnergy = ( Decay_Light.KE*Decay_Heavy.KE + Decay_Light.M*Decay_Heavy.KE + Decay_Heavy.M*Decay_Light.KE 
    -cos(Decay_Light.Theta - Decay_Heavy.Theta)
    *sqrt((Decay_Light.KE*Decay_Light.KE + 2*Decay_Light.M*Decay_Light.KE)
      *(Decay_Heavy.KE*Decay_Heavy.KE + 2*Decay_Heavy.M*Decay_Heavy.KE)) ) 
      /(Decay_Light.M + Decay_Heavy.M) + SeparationEnergy();
}

void Physics::CalcInvMassExcEnergy(){
  InvMassExcEnergy = Fragment.LV.M() - Fragment.M;
}


void Physics::CalcRecoilQValue(const double& dz /*mm*/, const double& tof /*ns*/){
  Ejectile.KE = 0.5*(Ejectile.M/90000)*(dz*dz)/(tof*tof);
  
  double t1 = (Beam.M/Fragment.M)*Beam.KE;
  double t2 = (Ejectile.M/Fragment.M)*Ejectile.KE;
  double t3 = (2/Fragment.M)*sqrt(Beam.KE*Ejectile.KE*Ejectile.M*Beam.M);  
  double d_heavy_ke = t1 + t2 + t3;

  RecoilQValue = Ejectile.KE + d_heavy_ke - Beam.KE;
}








#endif
