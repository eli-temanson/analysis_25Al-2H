/*

RxnTarget.cpp
A basic target unit for use in the SPANCRedux environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry 
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "RxnTarget.hpp"

RxnTarget::RxnTarget() {}
RxnTarget::~RxnTarget() {}

// /*Targets must be of known thickness*/
// RxnTarget::RxnTarget(double thick) {
//   thickness = thick;
// }

void RxnTarget::SetThickness(double thick){
  thickness = thick;
}

/*Set RxnTarget elements of given Z, A, S*/
void RxnTarget::SetElements(std::vector<int>& z, std::vector<int>& a, std::vector<int>& stoich) {
  Z = z;
  A = a;
  Stoich = stoich;
  
  eloss.SetTargetComponents(Z, A, Stoich);
}

/*Element verification*/
bool RxnTarget::ContainsElement(int z, int a) {
  for(unsigned int i=0; i<Z.size(); i++) {
    if( z == Z[i] && a == A[i]) return true;
  }
  return false;
}

/*Calculates energy loss for travelling all the way through the RxnTarget*/
double RxnTarget::getEnergyLossTotal(int zp, int ap, double startEnergy, double theta) {
  if(theta == PI/2.) return startEnergy;
  else if (theta > PI/2.) theta = PI - theta;
  return eloss.GetEnergyLoss(zp, ap, startEnergy, thickness/fabs(cos(theta)));
}

/*Calculates energy loss for travelling halfway through the RxnTarget*/
double RxnTarget::getEnergyLossHalf(int zp, int ap, double startEnergy, double theta) {
  if(theta == PI/2.) return startEnergy;
  else if (theta > PI/2.) theta = PI - theta;
  return eloss.GetEnergyLoss(zp, ap, startEnergy, thickness/(2.0*fabs(cos(theta))));
}

/*Calculates reverse energy loss for travelling all the way through the RxnTarget*/
double RxnTarget::getReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double theta) {
  if(theta == PI/2.) return finalEnergy;
  else if (theta > PI/2.) theta = PI - theta;
  return eloss.GetReverseEnergyLoss(zp, ap, finalEnergy, thickness/fabs(cos(theta)));
}

/*Calculates reverse energy loss for travelling half way through the RxnTarget*/
double RxnTarget::getReverseEnergyLossHalf(int zp, int ap, double finalEnergy, double theta) {
  if(theta == PI/2.) return finalEnergy;
  else if (theta > PI/2.) theta = PI - theta;
  return eloss.GetReverseEnergyLoss(zp, ap, finalEnergy, thickness/(2.0*fabs(cos(theta))));
}

/*Getter functions*/

double& RxnTarget::GetThickness() {
  return thickness;
}

int RxnTarget::GetNumberOfElements() {
  return Z.size();
}

int RxnTarget::GetElementZ(int index) {
  return Z[index];
}

int RxnTarget::GetElementA(int index) {
  return A[index];
}

int RxnTarget::GetElementStoich(int index) {
  return Stoich[index];
}
