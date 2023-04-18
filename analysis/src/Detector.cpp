#ifndef __Detector__cpp
#define __Detector__cpp

#include "Detector.hpp"

void Detector::Reset(){
  for(int i=0; i < Mult()+1; i++){
    Data[i].ERaw = Data[i].TRaw = Data[i].dERaw = ResetValue;
    Data[i].ChRaw = ResetValue;
    Data[i].E = ResetValue;
    Data[i].T = ResetValue;
  }
  _mult = 0; //reset the multiplicity
}

void Detector::Initialize(std::string s, const int& n){
  _name = s;
  _num = n;
  for(int i=0; i<NumOfCh(); i++){
    Calib[i].ChMap = i;
    Calib[i].A0 = Calib[i].E0 = Calib[i].T0 = Calib[i].TOffset = Calib[i].EOffset = 0.0;
    Calib[i].A1 = Calib[i].E1 = Calib[i].T1 = Calib[i].EScale = 1.0;

    Data[i].ERaw = Data[i].TRaw = Data[i].dERaw = ResetValue;
    Data[i].ChRaw = ResetValue;
  }
}


double Detector::ERaw(int i) const {
  return Data[i].ERaw;
}
double Detector::dERaw(int i) const {
  return Data[i].dERaw;
}
double Detector::TRaw(int i) const {
  return Data[i].TRaw;
}
double Detector::ChRaw(int i) const { 
  return Data[i].ChRaw;
}
double Detector::Ch(int i) const { 
  return Calib[Data[i].ChRaw].ChMap;
}


double Detector::E_Offset(int i) const {
  return ERaw(i) < 0 ? ResetValue : ERaw(i)*Calib[Data[i].ChRaw].EScale + Calib[Data[i].ChRaw].EOffset;
}
double Detector::E_LocalGain(int i) const {
  return E_Offset(i) < 0 ? ResetValue : E_Offset(i)*Calib[Data[i].ChRaw].A1 + Calib[Data[i].ChRaw].A0;
}
double Detector::E_LocalCalib(int i) const {
  return E_LocalGain(i)<0 ? ResetValue : E_LocalGain(i)*Calib[Data[i].ChRaw].E1 + Calib[Data[i].ChRaw].E0;
}
double Detector::E(int i) const { 
  return E_LocalCalib(i)<0 ? ResetValue : E_LocalGain(i)*e_calib_p2*E_LocalGain(i)*e_calib_p2 + E_LocalGain(i)*_ELin + _EShift;
}

// double Detector::dE_LocalGain(int i) const {
//   return (dERaw(i)<0) ? ResetValue : dERaw(i)*Calib[i].A1 + Calib[i].A0;
// }
// double Detector::dE_LocalCalib(int i) const {
//   return (dE_LocalGain(i)<0) ? ResetValue : dE_LocalGain(i)*Calib[i].E1 + Calib[i].E0;
// }
// double Detector::dE(int i) const { 
//   return (dE_LocalCalib(i)<0) ? ResetValue : dE_LocalCalib(i)*_ELin + _EShift;
// }

double Detector::T_Offset(int i) const {
  return TRaw(i) + Calib[Data[i].ChRaw].TOffset;
}
double Detector::T_LocalCalib(int i) const {
  return T_Offset(i)*Calib[Data[i].ChRaw].T1 + Calib[Data[i].ChRaw].T0;
}
double Detector::T(int i) const { 
  return T_LocalCalib(i)*_TLin + _TShift;
}



void Detector::SetCalibrations(VariableMap& VarMap){
  for(int i=0; i < NumOfCh(); i++){
    VarMap.GetParam(GetName()+".ch["+to_string(i)+"]",Calib[i].ChMap);
    VarMap.GetParam(GetName()+".a0["+to_string(i)+"]",Calib[i].A0);
    VarMap.GetParam(GetName()+".a1["+to_string(i)+"]",Calib[i].A1);
    VarMap.GetParam(GetName()+".e0["+to_string(i)+"]",Calib[i].E0);
    VarMap.GetParam(GetName()+".e1["+to_string(i)+"]",Calib[i].E1);
    VarMap.GetParam(GetName()+".t0["+to_string(i)+"]",Calib[i].T0);
    VarMap.GetParam(GetName()+".t1["+to_string(i)+"]",Calib[i].T1);
    VarMap.GetParam(GetName()+".e_scale["+to_string(i)+"]",Calib[i].EScale);
    VarMap.GetParam(GetName()+".e_offset["+to_string(i)+"]",Calib[i].EOffset);
    VarMap.GetParam(GetName()+".t_offset["+to_string(i)+"]",Calib[i].TOffset); 
  } 

  VarMap.GetParam(GetName()+"e_calib_p2",e_calib_p2);
  VarMap.GetParam(GetName()+".elowlimit",_LowELimit);
  VarMap.GetParam(GetName()+".ehighlimit",_HighELimit);
  VarMap.GetParam(GetName()+".elin", _ELin);
  VarMap.GetParam(GetName()+".eshift", _EShift );
  VarMap.GetParam(GetName()+".tlin",_TLin);
  VarMap.GetParam(GetName()+".tshift",_TShift);
}




void Detector::Print(){
  std::cout<<"\n--------------- Detector information for " << GetName() << "---------------\n";
  std::cout<<"Energy threshold: "<< _LowELimit <<" >< "<< _HighELimit <<  std::endl;
  std::cout<<"\nInput Data Structure \n";
  std::cout<<"ERaw    TRaw    ChRaw \n";
  for(int i=0; i<Mult(); i++){
    std::cout << Data[i].ERaw << "    " << Data[i].TRaw << "    " <<  Data[i].ChRaw << std::endl;
  }
  std::cout<<"Calibration for highest energy ch \n";
  std::cout<<"A0 = "<< Calib[Data[0].ChRaw].A0 <<"  A1 = "<< Calib[Data[0].ChRaw].A1 <<std::endl;
  std::cout<<"E0 = "<< Calib[Data[0].ChRaw].E0 <<"  E1 = "<< Calib[Data[0].ChRaw].E1 <<std::endl;
}




void Detector::InsertHit(const double& e, const double& t, const double& ch){
  if(e > _LowELimit && e < _HighELimit){
    Data[_mult] = {e,t,0.,ch};
    Data[_mult].E = E(_mult);
    Data[_mult].T = T(_mult);
    _mult += 1;
  }
}
void Detector::InsertHit(const double& e,const double& de, const double& t, const double& ch){
  if(e > _LowELimit && e < _HighELimit){
    Data[_mult] = {e,t,de,ch};
    Data[_mult].E = E(_mult);
    Data[_mult].T = T(_mult);
    _mult += 1;
  }
}
void Detector::SortByEnergy()
{
  // std::sort( std::begin(Data), std::end(Data), [](const Set& a, const Set& b) {return a.ERaw < b.ERaw;} );
  // Sort based on gain-matched and calibrated energy!
  // for(int i=0; i < Mult()+1; i++){
  //   Data[i].E = E(i); 
  // }
  // std::sort( std::begin(Data), std::begin(Data)+Mult(), [](const Set& a, const Set& b){ return a.ERaw > b.ERaw;} );
  
  std::sort( std::begin(Data), std::begin(Data)+Mult(), [](const Set& a, const Set& b){ return a.E > b.E;} );
}

void Detector::SortByTime()
{
  std::sort( std::begin(Data), std::begin(Data)+Mult(), [](const Set& a, const Set& b){ return a.T > b.T;} );
}



#endif
