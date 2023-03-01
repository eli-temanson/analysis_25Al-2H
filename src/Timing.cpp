#ifndef __Timing__cpp
#define __Timing__cpp

#include "Timing.hpp"

Timing::Timing(){};
Timing::~Timing(){};

void Timing::init(std::string name){
  _Name = name;
}

void Timing::Reset(){
  _T = ResetValue;
}

double Timing::TRaw() const{
  return _T;
}
double Timing::T_Offset() const{
  return TRaw() + _TOffset;
}
double Timing::T() const{ 
  return T_Offset()*_TLin + _TShift;
}
double Timing::TMod() const{
  return ( _T < 0 ? -200 : T() - RFperiod*std::floor(T()/RFperiod) ); 
  // return ( _T < 0 ? -200 : fmod( T(), RFperiod)); 
}
double Timing::TMod2() const{ //Modulus with an offset tzero
  return ( _T < 0 ? -200 : std::fmod( T()-_Tzero, RFperiod)); 
}
double Timing::TRel(double t) const{
  return ( _T < 0 ? -200 : std::fmod(std::abs(t-T()), RFperiod));
}


void Timing::SetCalibrations(VariableMap& VarMap){
  VarMap.GetParam(GetName()+".t0",_T0);
  VarMap.GetParam(GetName()+".t1",_T1);
  VarMap.GetParam(GetName()+".tlin",_TLin);
  VarMap.GetParam(GetName()+".tshift",_TShift);
  VarMap.GetParam(GetName()+".tzero",_Tzero);
  VarMap.GetParam(GetName()+".t_offset",_TOffset);
}


/* Insert Hit */
void Timing::InsertHit(const double& t){
  _T = t;
}


#endif
