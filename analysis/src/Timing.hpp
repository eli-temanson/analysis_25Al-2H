#ifndef __Timing__hpp
#define __Timing__hpp

#include "VariableMap.hpp"

#include <cmath>


class Timing{
protected:
  std::string _Name;
  double _T = 0;

  double ResetValue = -500.;
 
  double _TOffset=0.0;
  double _T0=0.0;
  double _T1=1.0;
  double _TLin=0.293;
  double _TShift=0.0;
  double _Tzero=0.0;
  
public:
  Timing(); //Constructor
  ~Timing(); //Destructor

  void init(std::string);
  void Reset();
  
  double RFperiod = 82.474227;
  double TRaw() const;
  double T_Offset() const;
  double T() const;
  double TMod() const;
  double TMod2() const;
  double TRel(double) const;

  std::string GetName() const{return _Name;}
  
  void InsertHit(const double&);  
  void SetCalibrations(VariableMap& VarMap);

};







#endif
