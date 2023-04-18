#ifndef __Global__hpp
#define __Global__hpp

#include "TROOT.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCutG.h"
#include "THashTable.h"

#include <iostream>
#include <string>

using namespace std;


/*2D histogram fill wrapper*/
inline void MyFill(THashTable *rootObj,string name,int binsx,double minx,double maxx,double valuex,
                                            int binsy,double miny,double maxy,double valuey){
  TH2F *histo = (TH2F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) {
    histo->Fill(valuex, valuey);
  } else {
    TH2F *h = new TH2F(name.c_str(), name.c_str(), binsx, minx, maxx, binsy, miny, maxy);
    h->Fill(valuex, valuey);
    rootObj->Add(h);
  }
}
/*1D histogram fill wrapper*/
inline void MyFill(THashTable *rootObj,string name,int binsx,double minx,double maxx,double valuex){
  TH1F *histo = (TH1F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) {
    histo->Fill(valuex);
  } else {
    TH1F *h = new TH1F(name.c_str(), name.c_str(), binsx, minx, maxx);
    h->Fill(valuex);
    rootObj->Add(h);
  }
}


inline bool InsideGate(TCutG* g,double x,double y){
  return g->IsInside(x,y);
}
inline bool InsideGate(double val,double low,double high){
  if(val > low && val < high){ return true; }
  else{ return false; }
}





#endif
