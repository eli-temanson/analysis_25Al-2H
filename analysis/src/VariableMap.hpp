/***********************************************************
 Sean Kuvin 2013
 RNVariableMap class

 This is the default parameter map that can be added to 
 any number of times depending on the implementation.  
 (i.e. load a silicon calibration file and the silicon will
 be automatically calibrated by the std::strings that that particular
 detector knows how to get such as "si_a.elin").  If no parameter
 is at some point loaded into the parameter map then the default values
 should be used.
***********************************************************/

#ifndef __VariableMap__hpp
#define __VariableMap__hpp

//C and C++ libraries.
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <memory>
#include <iterator>

class VariableMap{

protected:
  std::unordered_map<std::string,double> cal_map;
  std::unordered_map<std::string,double>::iterator iter;

public:
  VariableMap(){};
  ~VariableMap(){};
 

  void LoadParams(std::string filename){
    std::ifstream cal_file;
    std::string key;

    cal_file.open( filename.c_str() );

    if (!cal_file.is_open()){
      std::cout << "Could not open " << filename << std::endl;
      exit(0);
      return;
    }
    std::string temp_string;
    double value;
    do{
      temp_string="";
      value=0.0;
      cal_file>>temp_string>>value; 
      cal_map.insert(std::pair<std::string,double>(temp_string,value));

    }while(!cal_file.eof());
  }


  int GetParam(std::string param,float& var){
    iter=cal_map.begin();
    iter=cal_map.find(param);
    if(iter!=cal_map.end()){
      var=(float)iter->second;
      return 1;
    }
    else return 0;
  }

  int GetParam(std::string param,int& var){
    iter=cal_map.begin();
    iter=cal_map.find(param);
    if (iter!=cal_map.end()){
      var=(int)iter->second;
      return 1;
    }
    else return 0;
  }


  int GetParam(std::string param, double& var){
    iter=cal_map.begin();
    iter=cal_map.find(param);
    if (iter!=cal_map.end()){
      var=iter->second;
      return 1;
    }
    else return 0;
  }

  int GetParam(std::string param,bool& var){
    iter=cal_map.begin();
    iter=cal_map.find(param);
    if(iter!=cal_map.end()){
      var=(bool)iter->second;
      return 1;
    }
    else return 0;
  }


  void Print(){
    iter=cal_map.begin();
    while(iter!=cal_map.end()){
      std::cout<< iter->first <<" "<< iter->second <<"\n";
      iter++;
    }
  }

  int AddParam(std::string param,double var){
    cal_map.insert(std::pair<std::string,double>(param,var));
    return 1;
  }

  void ClearParams(){
    cal_map.clear();
  }

  
};


#endif
