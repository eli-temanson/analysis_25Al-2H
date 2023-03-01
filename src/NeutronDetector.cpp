#ifndef __NeutronDetector__cpp
#define __NeutronDetector__cpp

#include "NeutronDetector.hpp"

#define M_N 939.5653451
#define M_P 938.2720137 
#define M_C 11177.862342
#define M_AL 25133.14404
#define M_CU 59193.099

NeutronDetector::NeutronDetector(){};
NeutronDetector::~NeutronDetector(){};

void NeutronDetector::init(std::string name,const int& num){
  _Name=name;
  _Num=num;

  det.Initialize("neut.det",num);

  _T0=0.0;
  _T1=1.0;
  _TLin=1.0;
  _TShift=0.0;
  _ShortOffset=0.0;
  _APos=1;
  _HitPos.SetXYZ(0,0,0);
  _PhePol2=0.0;
  _PheLin=1.0;
  _PheShift=0.0;
  _ZRot=0.0;
  // _Radius=34.;
  // _Thickness=25.4;
  // _Threshold=0.05;
  // _TimingResolution=0.0;
  _QLong=0.0;
  _QShort=0.0;
  _PSD=0.0;
  _T=0.0;
  _TRel=0.0;
  //_PosGrid(Name()+".posgrid";num);
								  
  _Pos.SetXYZ(0,0,-228.7); //set default zpos
}


void NeutronDetector::Reset(){
  PosGrid.Reset();
  det.Reset();

  _QLong=0;
  _QShort=0;
  _TRel=0;
  _T=0;
  _PSD=0;
  
  // _Dt=0;
  // _ESum=0;
  // _TSim=0;
  // _Counter=0;
  // _CounterCarbon=0;
  // _ELost=0;
  _HitPos.SetXYZ(0,0,0);
}


void NeutronDetector::SetCalibrations(VariableMap& VarMap){
  det.SetCalibrations(VarMap);
  VarMap.GetParam(Name()+".tlin",_TLin);
  VarMap.GetParam(Name()+".tshift",_TShift);
  VarMap.GetParam(Name()+".t0",_T0);
  VarMap.GetParam(Name()+".t1",_T1);
  VarMap.GetParam(Name()+".shortoffset",_ShortOffset);
  
  double temp;
  if(VarMap.GetParam("neutarray.zpos",temp))
    _Pos.SetZ(temp);
  if(VarMap.GetParam(Name()+"%s.apos",temp)){
    _APos=(int)temp;
  }
  
  // VarMap.GetParam("neutarray.rotz",_ZRot);
  // VarMap.GetParam("all.radius",fRadius);
  // VarMap.GetParam("all.thickness",fThickness);
  // VarMap.GetParam("all.threshold",fThreshold);
  // VarMap.GetParam("all.timingresolution",fTimingResolution);
  // VarMap.GetParam(Form("%s.radius",Name()),fRadius);
  // VarMap.GetParam(Form("%s.thickness",Name()),fThickness);
  // VarMap.GetParam(Form("%s.threshold",Name()),fThreshold);
  // VarMap.GetParam(Form("%s.timingresolution",Name()),fTimingResolution);
  // VarMap.GetParam(Form("%s.phepol2",Name()),fPhePol2);
  // VarMap.GetParam(Form("%s.phelin",Name()),fPheLin);
  // VarMap.GetParam(Form("%s.pheshift",Name()),fPheShift);
  // _PosGrid.SetCalibrations(VarMap);

  // //set updated detector position in RNArray
  // PositionMap(_APos,_Pos,_Rotz);
}





//The IsANeutron/Gamma() methods show where I was hoping to go with the gating framework but I have not started using these yet
// int NeutronDetector::IsANeutron(){
//   TCutG *neut = (TCutG*)gROOT->GetListOfSpecials()->FindObject(Form("%s_neuts",Name()));
//   if(neut&& neut->IsInside(_QLong,_QShort)){
//     return 1;
//   }
//   else return 0;
// }
// int NeutronDetector::HitID(){
//   if(IsANeutron())
//     return 2;
//   else if(IsAGamma())
//     return 1;
//   else
//     return 0;
// }
// int NeutronDetector::IsAGamma(){
//   TCutG *gam = (TCutG*)gROOT->GetListOfSpecials()->FindObject(Form("%s_gammas",Name()));
//   if(gam && gam->IsInside(fQLong,fQShort))
//     return 1;
//   else return 0;
// }


// double NeutronDetector::T() const{
//   return fT>0 ? ((TLocal() * fTLin) + fTShift) : 0; //from TDC
// }

// double NeutronDetector::PosX() const{
//   double side=0.0;
//   double pos=0.0;
//   double posE=0.0;
//   for(unsigned int i=0;i<PosGrid.Mult();i++){
//     if(PosGrid.Ch(i)==0 || PosGrid.Ch(i)==1){
//       side = -1;
//     }
//     else if(PosGrid.Ch(i)==2 || PosGrid.Ch(i)==3){
//       side = 1;
//     }
//     pos += (side *PosGrid.E(i));
//     posE += PosGrid.E(i);
//   }
//   return pos/posE;
// }

// double NeutronDetector::PosY() const{
//   double side=0.0;
//   double pos=0.0;
//   double posE=0.0;
//   for(unsigned int i=0;i<fPosGrid.Mult();i++){
//     if(fPosGrid.Ch(i)==0 || fPosGrid.Ch(i)==3){
//       side = -1;
//     }
//     else if(fPosGrid.Ch(i)==1 || fPosGrid.Ch(i)==2){
//       side = 1;
//     }
//     pos += (side * fPosGrid.E(i));
//     posE += fPosGrid.E(i);
//   }
//   return pos/posE;
// }

// double NeutronDetector::SumPosE()const{
//   double e=0.0;
//   for(unsigned int i=0;i<fPosGrid.Mult();i++){
//     e += fPosGrid.E(i);
//   }
//   return e;
// }






//just use Z distance (+1/2 detector thickness) to approximate neutron as on axis
// double NeutronDetector::nKE(double tof) const{
//   return ((.5*939.565*ZEst()*ZEst())/(90000*tof*tof));
// }

// //use actual path of neutron to get energy
// double NeutronDetector::nKE_R(double tof) const{
//   double r = ZEst()/TMath::Cos(fPos.Theta());
//   return ((.5*939.565*r*r)/(90000*tof*tof));
// }


// double NeutronDetector::Q_value_est(double tof,double mass1,double mass2,
// 				                            double beam_e,double& hi_KE,double& Q_val){
//   double ne = nKE(tof);
//   hi_KE = ((mass1*beam_e/mass2)+(939.565*ne/mass2)+((2/mass2)*sqrt(beam_e*ne*mass1*939.565)));
//   Q_val = ne+hi_KE-beam_e;
  
//   return Q_val;
// }



// void NeutronDetector::InsertPSDHit(const double& q_long,const double& q_short,const double& t){
//   _QLong=q_long;
//   _QShort=q_short;
//   _T=t;
  
//   if(_QShort>0 && _QLong>0) 
//     _PSD = QShortOffset()/_QLong;
// }


// double NeutronDetector::CalculateTRel(const double &tfirst){
//   if(_T)
//     _TRel = T()-_first;
//   return _TRel;
// }


// double NeutronDetector::QkeVee() const {
//   return (_QLong*_QLong*_PhePol2 + _QLong*_PheLin + _PheShift);
// }




//-------------------------------------------------------------------------------
//
//
// Neutron detector array begins here
//
//
//-------------------------------------------------------------------------------
// NeutronDetectorArray::NeutronDetectorArray(const TString & name):BaseClass(name),
//   fT_first=0.0,
//   fDetfirst(-1),
//   fT_mult=0.0,
//   fMult=0.0,
//   fPos(16,TVector3(0,0,0)),
//   fQ_long(16,0),
//   fPSD(16,0),
//   fT(16,0),
//   fDetlist(16,-1)
// {  
// }



// int NeutronDetectorArray::ReconstructHits(NeutCollection& in){
  
//   int cref=0;
//   for(NeutCollectionRef it=in.begin();it!=in.end();it++){  
//     if((*it).QLong()>0){
//       double ql=(*it).QLong();
//       double qs=(*it).QShortOffset();
//       TVector3 pos=(*it).GetPosVect();
//       double t=(*it).T();
//       InsertHit(ql,qs,t,pos,cref);
//     }
//     cref++;
//   }

//   return 1;
// }


// void NeutronDetectorArray::Reset(){
//   for(int i=0;i<fMult;i++){
//     fPos[i].SetXYZ(0,0,0);
//     fQ_long[i]=0;
//     fPSD[i]=0;
//     fDetlist[i]=-1;
//     fT[i]=0;
//   }
//   fMult=0;
//   fT_mult=0;
//   fT_first=0;
//   fDetfirst = -1;

//   return;
// }


// int NeutronDetectorArray::InsertHit(const double& q_long,const double& q_short,const double& q_t,
//                                     const TVector3& pos,const int& index){
  
//   if (!q_long>0){
//     return -1;
//   }

//   int i,j;
//   double psd = (q_short/q_long);

//   /* sorted by energy */
//   for (i=(int)fMult-1;i>=0;i--){
//     if(q_long < fQ_long[i]){
//       break;
//     }
//   }
  
//   /* element i+1 is at the position for ch so make room for it */
//   for (j=(int)_Mult-1; j>i; j--){
//     _Detlist[j+1]=_Detlist[j];
//     _Qlong[j+1]=_Qlong[j];
//     _PSD[j+1]=_PSD[j];
//     _T[j+1]=_T[j];
//     _Pos[j+1]=_Pos[j];
//   }


//   // and shove it in
//   _Detlist[i+1]=index;
//   _Qlong[i+1]=q_long;
//   _PSD[i+1]=psd;
//   _Pos[i+1]=pos;
//   _T[i+1]=q_t;
//   _Mult += 1;

//   return (i+1);
// }





// bool NeutronDetector::inDet(const TVector3& v, const TVector3& beamspot){
//   /* first see if it intersects the plane of the detector */
//   TVector3 normv(0,0,_Pos.Z());
//   double udotnormv = (v.Unit()).Dot(normv);
//   if(udotnormv <= 0.){return false;}

//   /*Find distance from target to interesection point of detector plane 
//     use line-plane intersection formula */
//   TVector3 posvect = GetPosVect();
//   //  double vdist = (fPos.Unit()).Dot(posvect)/udotnormv;
//   double vdist = (normv).Dot(posvect)/udotnormv;
//   if(vdist <= 0.){return false;}
  
//   //vector from target to interesection point of detector plane with
//   //magnitude equal to distance
//   TVector3 v_to_det(v);
//   v_to_det.SetMag(vdist);

//   //create vector from detector origin to interaction point on plane of detector
//   TVector3 ch_vect = v_to_det + beamspot - posvect;

//   //see if it is within the detector area
//   if((ch_vect.Mag()) > _Radius){
//     return false;
//   }
//   return true;
// }




// int NeutronDetector::NeutIn(TLorentzVector nLV,double& t,double& e,const TVector3& beamspot){
//   TLorentzVector inLV = nLV;
//   double pz=nLV.Pz(), px=nLV.Px(),py=nLV.Py();
//   double psqr=px*px+py*py+pz*pz;
//   double tof=(fPos.Z())*M_N/(pz*3*100);//tof to front of detector
//   double x_pos=(px*tof*300/(M_N))-fPos.X() + beamspot.X();
//   double y_pos=(py*tof*300/(M_N))-fPos.Y() + beamspot.Y();
//   double radial_pos=sqrt(x_pos*x_pos+y_pos*y_pos);  
//   double dx=0.0,dy=0.0,dz=0.0;
//   double nKE=0.0;
//   double macrosigma=0.0;
//   double step=0.0;

//   fHitPos.SetXYZ(x_pos+fPos.X(),y_pos+fPos.Y(),fPos.Z());
//   stepcounter = 0;
//   steptime = 0; 
//   z_pos = fThickness;
//   fDt = 0; //time to first reaction

//   do{
//     nKE=inLV.E()-inLV.M(); //get current neutron incident energy
//     if(nKE<fThreshold){ //check if it is above threshold
//       break;
//     }
//     px = inLV.Px();
//     py = inLV.Py();
//     pz = inLV.Pz();
    
//     //Take a step following an exponential distribution
//     macrosigma = cross_sections::neutron_pterphenyl_total_macro(nKE);
//     step = pdf::exponential(macrosigma);
//     steptime += fabs(step*M_N/(inLV.P()*3*100));
//     stepcounter++;

//     //make sure step is still within the detector dimensions
//     dx=step*px/(sqrt(psqr));
//     dy=step*py/(sqrt(psqr));
//     dz=step*pz/(sqrt(psqr));
    
//     x_pos += dx;
//     y_pos += dy;
//     z_pos += dz;//move to start of next sector.
//     radial_pos = sqrt(x_pos*x_pos+y_pos*y_pos);
//     if((radial_pos > fRadius || z_pos < 0 || z_pos > fThickness)){
//       break;
//     }

//     //we already found interaction point, now we need to decide what type of interaction took place
//     if(global::myRnd.Rndm()<cross_sections::neutron_pterphenyl_hydrogen_macro(nKE) / cross_sections::neutron_pterphenyl_total_macro(nKE)){
//       H_hit(inLV);
//     }
//     else{
//       C_hit(inLV);
//     }
//   }while (radial_pos <= fRadius && z_pos >= 0 && z_pos <= fThickness); //redundant at the moment
  
//   if(fESum>fThreshold){
//     e = fESum;
//     t = tof+ fDt;
//     if(fTimingResolution>0){
//       fTSim = t - global::myRnd.Gaus(0,fTimingResolution / 2.355); //fTimingResolution in FWHM to sigma
//     }
//     else{
//       fTSim=t;
//     }
//     return true;
//   }
//   else{
//     e=0;
//     t=0;
//     return false;
//   }
  
// }	  


// double NeutronDetector::H_hit(TLorentzVector& inLV){
  
//   //Create the Neut+Target LV
//   TLorentzVector neut_LVcopy(0,0,inLV.P(),inLV.E());
//   TLorentzVector Target(0.,0.,0.,M_P);
//   TLorentzVector Before(Target + neut_LVcopy);

//   //Boost into CM frame
//   TVector3 boostv = (-1)*Before.BoostVector();
//   neut_LVcopy.Boost(boostv); 

//   //random CM_Angles determines kinematics
//   double theta = acos(-1. + 2.*global::myRnd.Rndm());
//   double phi = 2.* 3.14 * global::myRnd.Rndm(); //isotropic CM
//   neut_LVcopy.SetTheta(theta);
//   neut_LVcopy.SetPhi(phi);


//   //Set CM LV and boost back to lab frame
//   neut_LVcopy.Boost(-1*boostv);//boost back to lab frame


//   //Rotate vector back from z-axis
//   neut_LVcopy.RotateY(inLV.Theta());
//   neut_LVcopy.RotateZ(inLV.Phi());

  
//   //Get Energy Loss
//   double nEdep = inLV.E()-neut_LVcopy.E();
//   if(nEdep>fThreshold){
//     fESum += nEdep;
//     //increment hit counter and if its the first hit, save the steptime
//     fCounter++;
//     if(fCounter==1){
//       fDt=steptime;
//     }
//   }
//   else{
//     fELost += nEdep;
//   }
//   //return the new E/angle to the original neut Lorentz Vector
//   inLV = neut_LVcopy;

//   return theta;
  
// }

// double NeutronDetector::C_hit(TLorentzVector& inLV){
  
//   //Create the Neut+Target LV
//   TLorentzVector neut_LVcopy(0,0,inLV.P(),inLV.E());
//   TLorentzVector Target(0.,0.,0.,M_C);
//   TLorentzVector Before(Target + neut_LVcopy);

//   //Boost into CM frame
//   TVector3 boostv = (-1)*Before.BoostVector();
//   neut_LVcopy.Boost(boostv);

//   //increment Carbon hit counter to see how much energy gets lost
//   fCounterCarbon++;

//   //Calculate En CM + random CM_Angles
//   double theta = acos(-1. + 2.*global::myRnd.Rndm());
//   double phi = 2.* 3.14 * global::myRnd.Rndm(); //isotropic CM

//   //Set CM LV and boost back to lab frame
//   neut_LVcopy.SetTheta(theta);
//   neut_LVcopy.SetPhi(phi);
//   neut_LVcopy.Boost(-1*boostv);//boost back to lab frame


//   //Rotate vector back from z-axis
//   neut_LVcopy.RotateY(inLV.Theta());
//   neut_LVcopy.RotateZ(inLV.Phi());

    
//   //Get Energy Loss
//   fELost += inLV.E()- neut_LVcopy.E();//add what neut loses

//   //return the new E/angle to the original neut Lorentz Vector
//   inLV = neut_LVcopy;

//   return theta;  
// }



// int NeutronDetector::NeutInAfterScattering(TLorentzVector &nLV,double& t, double & e){
//   TLorentzVector inLV = nLV;
//   double dtnip=0.0,step=0.0;
//   const double * cloc = gGeoManager->GetCurrentPoint();
//   fHitPos.SetXYZ(cloc[0],cloc[1],cloc[2]);
//   const double * cdir = gGeoManager->GetCurrentDirection();
//   int detnode = gGeoManager->FindNode()->GetNumber(); //get the current node
//   int currentnode = gGeoManager->FindNode()->GetNumber(); //get the current node
  
//   //which corresponds to this detector
//   steptime = 0; //resetting step timer
  
//   do{
//     double nKE=inLV.E()-inLV.M(); //get current neutron incident energy
    
//     //Take a step following an exponential distribution
    
//     double macrosigma = cross_sections::neutron_pterphenyl_total_macro(nKE);
//     dtnip = pdf::exponential(macrosigma); //to cm
//     gGeoManager->SetCurrentDirection(inLV.Px(),inLV.Py(),inLV.Pz());
//     const double *cdir2 = gGeoManager->GetCurrentDirection();
//     const double *cloc2 = gGeoManager->GetCurrentPoint();
//     gGeoManager->FindNode();
//     gGeoManager->FindNextBoundary();
//     step = gGeoManager->GetStep();
    
//     //make sure step is still within the detector dimensions
//     if(nKE<fThreshold){ //check if it is above threshold
//       gGeoManager->SetStep(step+(0.01/inLV.P()));
//       gGeoManager->Step(kFALSE);
//       steptime += fabs(step*M_N/300); // dR = step * P()
//       stepcounter++; //if not, push it out of the detector and break
//       break;
//     }
    
//     if(dtnip>step*inLV.P()){
//       gGeoManager->SetStep(step+(0.01/inLV.P()));
//       gGeoManager->Step(kFALSE);
//       steptime += fabs(step*M_N/300); // dR = step * P()
//       stepcounter++;
//     }else{
//       //we already found interaction point, now we need to decide what type of interaction took place
//       steptime += fabs(dtnip*M_N/(inLV.P()*300));
//       stepcounter++;
//       gGeoManager->SetStep(dtnip/inLV.P());
//       gGeoManager->Step();
      
//       if(global::myRnd.Rndm()<cross_sections::neutron_pterphenyl_hydrogen_macro(nKE) /    cross_sections::neutron_pterphenyl_total_macro(nKE)){
	    
//         H_hit(inLV);
      
//       }else{
// 	      C_hit(inLV);
//       }
//     }
    
//     if(gGeoManager->IsOutside()){
//       break;
//     }

//     currentnode = gGeoManager->FindNode()->GetNumber(); //get the current node
//   }while(currentnode==detnode); 
  

//   if(fESum>fThreshold){
//     e = fESum;
//     t += fDt;
//     if (fTimingResolution>0){
//       fTSim = t - global::myRnd.Gaus(0,fTimingResolution / 2.355); //fTimingResolution in FWHM to sigma
//     }
//     else{
//       fTSim=t;
//     }
//     return true;
//   }
//   else{
//     e=0;
//     t=0;
//     return false;
//   }
  
// }	  



// int PositionMap(int slot,TVector3 & pos, double rot){
//   double z=pos.Z();
  
//   if(slot==1)pos.SetXYZ(38.1,-152.4,z);
//   else if(slot==2)pos.SetXYZ(-38.1,-152.4,z);
//   else if(slot==3)pos.SetXYZ(114.3,-76.2,z);
//   else if(slot==4)pos.SetXYZ(38.1,-76.2,z);
//   else if(slot==5)pos.SetXYZ(-38.1,-76.2,z);
//   else if(slot==6)pos.SetXYZ(-114.3,-76.2,z);
//   else if(slot==7)pos.SetXYZ(150.622,0,z);
//   else if(slot==8)pos.SetXYZ(74.422,0,z);
//   else if(slot==9)pos.SetXYZ(-74.422,0,z);
//   else if(slot==10)pos.SetXYZ(-150.622,0,z);
//   else if(slot==11)pos.SetXYZ(114.3,76.2,z);
//   else if(slot==12)pos.SetXYZ(38.1,76.2,z);
//   else if(slot==13)pos.SetXYZ(-38.1,76.2,z);
//   else if(slot==14)pos.SetXYZ(-114.3,76.2,z);
//   else if(slot==15)pos.SetXYZ(38.1,152.4,z);
//   else if(slot==16)pos.SetXYZ(-38.1,152.4,z);
//   else return 0;
  
  
//   pos.RotateZ(rot*TMath::DegToRad());
//   return 1;
  
// }








//-------------------------------------------------------------------------------------
//
// Geometry Package Methods                       
//
//-------------------------------------------------------------------------------------

// TGeoVolume* NeutronDetector::CreateCrystal(const int& node,const TVector3 & shift){
//   if(RNArray::gRNVolume){
//    return CreateCrystal(RNArray::gRNVolume,node,shift);
//   }   
//   else{
//     std::cout<<"use RNArray::CreateRNVolume() then NeutronDetector:CreateCrystal(node,shift)\n Or, use NeutronDetector::CreateCrystal(parent,node,shift)"<<std::endl;
//     return 0;
//   }
// }

// TGeoVolume* NeutronDetector::CreateCrystal(TGeoVolume*parent,const int& node,const TVector3 & shift){
//   TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
//   TGeoMedium *pterph = new TGeoMedium("Vacuum",1,mat);
//   TGeoVolume* crystal = gGeoManager->MakeTube(Name(),pterph,0,fRadius,fThickness/2);
//   parent->AddNode(crystal,node,new TGeoTranslation(fPos.X(),fPos.Y(), fPos.Z() - (fThickness/ 2.)));
//   return crystal;
// }

// TGeoVolume* NeutronDetector::CreateSilicone(const int& node,const TVector3 & shift){
//   if(RNArray::gRNVolume){
//    return CreateCrystal(RNArray::gRNVolume,node,shift);
//   }   
//   else{
//     std::cout<<"use RNArray::CreateRNVolume() then NeutronDetector:CreateCrystal(node,shift)\n Or, use NeutronDetector::CreateCrystal(parent,node,shift)"<<std::endl;
//     return 0;
//   }
// }


// TGeoVolume* NeutronDetector::CreateSilicone(TGeoVolume*parent,const int& node,const TVector3 & shift){
//   TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
//   TGeoMedium *silmed = new TGeoMedium("Vacuum",1,mat);
//   TGeoVolume* crystal = gGeoManager->MakeBox(Form("silicone_%s",Name()),silmed,32,32,4.75);
//   TGeoTranslation* trmove = new TGeoTranslation(fPos.X(),fPos.Y(),fPos.Z()-60);
//   TGeoRotation* rotmove = new TGeoRotation("silrotmove",0,0,45);
//   TGeoCombiTrans c1(*trmove,*rotmove);
//   TGeoHMatrix * hmsilicone = new TGeoHMatrix(c1);
//   parent->AddNode(crystal,node,hmsilicone);
//   return crystal;
// }




//-------------------------------------------------------------------------------------
//
// Previously RNArray namespace                        
//
//-------------------------------------------------------------------------------------

// TVector3 gNeutTrack(0,0,0); //To keep track of neutron position

// void ResetTracks(){
//   gNeutTrack.SetXYZ(0,0,0);
// }


// TGeoVolume *CreateRNVolume(){
//   TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
//   TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);
//   gRNVolume = gGeoManager->MakeTube("resoneut",med,0,210.0,600); //60 cm to track the entire array

//   return gRNVolume;  //for neutron scattering in the volume, we will make this the top volume . gGeoManager->SetTopVolume(gRNVolume);
// }



// TGeoVolume *AluminiumHousing::CreateORING(TGeoVolume* parent,const int&node,const TVector3&shift){
//   TGeoMaterial *mat = new TGeoMaterial("rubber",0,0,0);
//   TGeoMedium *med = new TGeoMedium("rubber",1,mat);
//   TGeoVolume* oring = gGeoManager->MakeTube("ORING",med,135.45,145.45,1.16); //60 cm to track the entire array
//   parent->AddNode(oring,node,new TGeoTranslation(fPos.X(),fPos.Y(),fPos.Z()+2.54+13+101.72+1.16));
//   oring->SetLineColor(kBlack);
//   return oring;
// }



// //I am implementing the Housing as if it was a detector itself, since I hope to keep track of neutrons which get deflected from the housing as opposed to
// //neutrons which are incident directly on the detector area.
// AluminiumHousing::AluminiumHousing(const std::string& name):fHitPos(TVector3(0,0,0)),
//                   fPos(TVector3(0,0,-228.7)),
//                   fAluminiumCounter(0),
//                   fELost(0),
//                   fOuterrad(200.32),
//                   fInnerrad(19.),
//                   fMinThreshold(0.01),
//                   fThickness(78.5),
//                   fThreshold(0.01),
//                   fRotz(0)
// {
//   char* pPath;
//   pPath = getenv ("ANALYSIS");
//   if(pPath!=NULL){
//     angdis[0].Init(Form("%s/input/AngularDistributions/NeutronScattering/Aluminium_400keV.dat",pPath));
//     angdis[1].Init(Form("%s/input/AngularDistributions/NeutronScattering/Aluminium_700keV.dat",pPath));
//     angdis[2].Init(Form("%s/input/AngularDistributions/NeutronScattering/Aluminium_1600keV.dat",pPath));
//     angdis[3].Init(Form("%s/input/AngularDistributions/NeutronScattering/Aluminium_3000keV.dat",pPath));
//   }
// }


// void AluminiumHousing::Reset(){
//   fAluminiumCounter = 0;
//   fAbsorberCounter = 0;
//   fELost =0;
//   fHitPos.SetXYZ(0,0,0);
//   fHitID = 0;
// }

// void AluminiumHousing::SetCalibrations(VariableMap& VarMap){
//   double temp;
//   if(VarMap.GetParam("neutarray.zpos",temp))
//     fPos.SetZ(temp);
//   VarMap.GetParam("neutarray.rotz",fRotz);
  
  
//   //set updated detector position in RNArray
  
// }

// TGeoVolume * AluminiumHousing::CreateVolume(TGeoVolume*parent){
//   TVector3 shift(fPos.X(),fPos.Y(),fPos.Z() - (fThickness/2)); //center the volume 
//   return CreateVolume(parent,kHOUSING,shift);
// }

// TGeoVolume * AluminiumHousing::CreateVolume(TGeoVolume*parent ,const int& node,const TVector3&shift){
//   TVector3 pos(0,0,0);
//   TGeoTube* stube = new TGeoTube("alH",fInnerrad,fOuterrad,(fThickness/2));
//   double origin[3] = {0,0,32-fThickness/2};
//   TGeoBBox *sbox = new TGeoBBox("alB",34,34,32,origin);
//   TGeoTube *stube2 = new TGeoTube("alC",0,34,(fThickness-2.54)/2);
//   TGeoTranslation *trans[16];
//   TGeoTranslation *trans2[16];
//   TGeoTranslation *trans33 = new TGeoTranslation("transmain",0,0,2.54);
//   trans33->RegisterYourself();
//   for(int i=0; i<16;i++){
//     PositionMap(i+1,pos,0);//We don't rotate here, because the entire housing needs to be rotated together
//     trans[i] = new TGeoTranslation(Form("trans%d",i),pos.X(),pos.Y(),0);
//     trans2[i] = new TGeoTranslation(Form("trans%d_2",i),pos.X(),pos.Y(),2.54/2);
//     trans[i]->RegisterYourself();    
//     trans2[i]->RegisterYourself();
//   }

//     TGeoCompositeShape *compshape = new TGeoCompositeShape("compshape","(alH:transmain - (alB:trans1 + alB:trans2 + alB:trans3 + alB:trans4 + alB:trans5 + alB:trans6 + alB:trans7 + alB:trans8 + alB:trans9 + alB:trans10 + alB:trans11 + alB:trans12 + alB:trans13 + alB:trans14 + alB:trans15 + alB:trans0) - (alC:trans1_2 + alC:trans2_2 + alC:trans3_2 + alC:trans4_2 + alC:trans5_2 + alC:trans6_2 + alC:trans7_2 + alC:trans8_2 + alC:trans9_2 + alC:trans10_2 + alC:trans11_2 + alC:trans12_2 + alC:trans13_2 + alC:trans14_2 + alC:trans15_2 + alC:trans0_2))");

//     TGeoRotation *rot1 = new TGeoRotation("hrot",0,0,fRotz);
//     TGeoTranslation *tr1 = new TGeoTranslation("htrans",shift.X(),shift.Y(),shift.Z());



//     TGeoVolume *housing = new TGeoVolume("housing",compshape);
//     parent->AddNode(housing,node,new TGeoCombiTrans(*tr1,*rot1));
//     return housing; 
    

// }

// TGeoVolume* AluminiumHousing::CreateFrontVolume(TGeoVolume*parent){
//   TVector3 shift(0,0,fPos.Z() + 114.72 / 2.+ 2.54); 
//   return CreateFrontVolume(parent,kFHOUSING,shift);
  
// }

// TGeoVolume* AluminiumHousing::CreateFrontVolume(TGeoVolume*parent ,const int& node,const TVector3&shif){    
//   TGeoPgon* pgon = new TGeoPgon("pgon",0,360,8,2);
//   pgon->DefineSection(0,50.86,0,180.45);  
//   pgon->DefineSection(1,-50.86,0,180.45);  
//   TGeoVolume*pgvol = new TGeoVolume("pgvol",pgon);
//   TGeoTube * frontin = new TGeoTube("frontin",0,101.45+79-13.7,42.69);
//   TGeoTube * frontopen = new TGeoTube("frontopen",0,101.45,8.17);
//   TGeoTube * mid = new TGeoTube("mid",101.45+79-13.7,200.32,6.5);
//   TGeoTranslation*trans11 = new TGeoTranslation("t11",0,0,6.5);
//   TGeoTranslation*trans12 = new TGeoTranslation("t12",0,0,50.86+6.5-8.17);
//   TGeoTranslation*trans13 = new TGeoTranslation("t13",0,0,-50.86); 
//   TGeoTranslation*trans14 = new TGeoTranslation("t14",0,0,6.5-8.17); 
  
  
//   trans11->RegisterYourself();
//   trans12->RegisterYourself();
//   trans13->RegisterYourself();
//   trans14->RegisterYourself();
  
  
//   TGeoCompositeShape* compshape = new TGeoCompositeShape("compshape"," pgon:t11+ mid:t13 - frontin:t14 - frontopen:t12");      
  
//   TGeoVolume * fhousing = new TGeoVolume("fhousing",compshape);
//   parent->AddNode(fhousing,node,new TGeoTranslation("hTRANs",shift.X(),shift.Y(),shift.Z()));
//   return fhousing; 
  
// }

// TGeoVolume* AluminiumHousing::CreateTargetVolume(TGeoVolume*parent){
//   double targ_zpos = (fPos.Z() + 2.54 + 13 + 101.72 + 2.32 + 13 + 63.16);
//   return CreateTargetVolume(parent,kTHOUSING,TVector3(0,0,targ_zpos));
// }

// TGeoVolume* AluminiumHousing::CreateTargetVolume(TGeoVolume*parent ,const int& node,const TVector3&shift){
//   TGeoTube * targend = new TGeoTube("targend",0,144.8,6.5);
//   TGeoTranslation*targfront =  new TGeoTranslation("targfront",0,0,76.16-6.5);
//   TGeoTranslation*targback =  new TGeoTranslation("targback",0,0,-76.16+6.5);
//   TGeoBBox *targbox = new TGeoBBox("targbox",93.75,93.75,(152.32-26)/2);
//   TGeoBBox *targinside = new TGeoBBox("targinside",75.9,75.9,(152.32/2));
//   targfront->RegisterYourself();
//   targback->RegisterYourself();
  
//   TGeoCompositeShape* tcshape = new TGeoCompositeShape("tcshape","targbox+targend:targfront+targend:targback-targinside");      
//   TGeoVolume * thousing =new TGeoVolume("targethousing",tcshape);
//   parent->AddNode(thousing,node,new TGeoTranslation(shift.X(),shift.Y(),shift.Z())); 
  
//   return thousing; 
  
// }

// TGeoVolume * AluminiumHousing::CreateBackingVolume(TGeoVolume*parent){
//   return CreateBackingVolume(parent,kBHOUSING,TVector3(0,0,-78.5-12.805));
  
// }


// TGeoVolume * AluminiumHousing::CreateBackingVolume(TGeoVolume*parent ,const int& node,const TVector3&shift){
//   TGeoBBox * outbox =new TGeoBBox("outbox",38.1,38.1,12.805);
//   TGeoBBox * inbox =new TGeoBBox("inbox",34,34,11.26);
//   TGeoTranslation* inboxtrans = new TGeoTranslation("inboxtrans",0,0,12.805-11.26);
//   TGeoRotation* rotmove = new TGeoRotation("rotmove",0,0,45);
//   inboxtrans->RegisterYourself();
//   rotmove->RegisterYourself();
//   TGeoCombiTrans *c1 = new TGeoCombiTrans(*inboxtrans,*rotmove);
//   c1->SetName("boxmove");
//   c1->RegisterYourself();
//   TGeoCompositeShape * cshape = new TGeoCompositeShape("backbox","(outbox:rotmove-inbox:boxmove)");
//   TGeoVolume * backbox = new TGeoVolume("backbox",cshape);

//   for(int i=0;i<16;i++){
//     TVector3 pos(0,0,fPos.Z());
//     RNArray::PositionMap(i+1,pos,fRotz);
//     TGeoTranslation* trmove = new TGeoTranslation(pos.X()+shift.X(),pos.Y()+shift.Y(),pos.Z()+shift.Z());
//     parent->AddNode(backbox,node+i,trmove);
    
//   }
//   return backbox;
// }


// Double32_t AluminiumHousing::GetNextInteractionPoint(double nKE /*Lab in MeV*/){
//   double macrosigma = cross_sections::neutron_aluminium_macro(nKE);
//   double step = pdf::exponential(macrosigma);
//   return step;
// }

// Double32_t AluminiumHousing::GetAbsorberInteractionPoint(double nKE /*Lab in MeV*/){
//   double macrosigma = cross_sections::neutron_copper_macro(nKE);
//   double step = pdf::exponential(macrosigma);
//   return step;
// }




// int AluminiumHousing::NeutIn(TLorentzVector nLV,double& t,double& e,const TVector3& beamspot){
//   TLorentzVector inLV = nLV;
//   double pz=nLV.Pz(), px=nLV.Px(),py=nLV.Py();
//   double psqr=px*px+py*py+pz*pz;
//   double step(0),dtnip(0);
//   steptime = 0;
//   stepcounter = 0;
//   fgIsDirect = 0;
//   fgNeutHit = 0;
//   //We get the hit position on the front of the housing, but the zposition is chosen
//   //considering that the "World" is now the resoneut housing volume which consists of the 
//   //detectors and the aluminium frame. the origin of this world is the center of the resoneut housing.
//   //therefore the track Z starts at the front of the housing (fThickness/2) and we subtract .1 mm 
//   //to ensure that the tracks are inside of the frame.  
//   gNeutTrack.SetXYZ(0,0,0); //be aware that these values are in cm 
//   gGeoManager->InitTrack(gNeutTrack.X(),gNeutTrack.Y(),gNeutTrack.Z(),px,py,pz); //be aware that these values are in cm
//   gGeoManager->FindNode();
  
//   do{
//     double nKE=inLV.E()-inLV.M(); //get current neutron incident energy
//     if(nKE<fThreshold){break;} //check if it is above threshold
  
//     if(gGeoManager->IsOutside()){return 0;}//particle has left resoneut

//     int nodenum = gGeoManager->FindNode()->GetNumber();

//     if(nodenum == kHOUSING || nodenum == kFHOUSING || nodenum == kTHOUSING || (nodenum >= kBHOUSING && nodenum < kBHOUSING+16)){
// dtnip = GetNextInteractionPoint(nKE);	
// gGeoManager->FindNextBoundary();
// step = gGeoManager->GetStep();

// if(step*inLV.P()>dtnip){
//   gGeoManager->SetStep(dtnip/inLV.P()); //only step to next interaction point
//   gGeoManager->Step();
//   Al_hit(inLV);
  
//   if(fAluminiumCounter==1){
//     const double *cloc = gGeoManager->GetCurrentPoint();
//     fHitPos.SetXYZ(cloc[0],cloc[1],cloc[2]);
//     fHitID = nodenum; //store id for which part of housing was hit 
//   }
  
//   gGeoManager->SetCurrentDirection(inLV.Px(),inLV.Py(),inLV.Pz());
//   steptime += fabs(dtnip*M_N/(inLV.P()*300));
//   stepcounter++;
//   continue;
// }
// else{//step all the way to boundary and cross
//   gGeoManager->SetStep(step+(0.01/inLV.P())); //cross by 0.01 mm
//   gGeoManager->Step(kFALSE);
//   steptime += fabs(step*M_N/300); //dR = step *inLV.P()
//   stepcounter++;
//   continue;
// }
//     }
    
    
//     if(nodenum >= 17 && nodenum<NEUTNUM+17){
// int index = nodenum - 17;
// double savetime = steptime;
// t=steptime; //set time 

// if(RNROOT::neut[index].NeutInAfterScattering(inLV,t,e)){
//   fgNeutHit = 1;
//   steptime = t;
// }//if a neutron isnt detected the t returned by this method is 0
// else{
//   steptime += savetime; //steptime from the neutron detector counter //this is getting more hack-y
// }
// //return 1; break after finding new detector
// continue;
//     }

    
//     if(nodenum == 1){
// gGeoManager->FindNextBoundary(0.5); //take small steps
// step = gGeoManager->GetStep();
// gGeoManager->SetStep(step+(0.01/inLV.P())); //cross by 0.01 mm
// gGeoManager->Step(kFALSE);
// steptime += fabs(step*M_N/300); //dR = step *inLV.P()
// stepcounter++;
// continue;
//     }

//     if(nodenum >= kSILICONE && nodenum < (NEUTNUM + kSILICONE)){
// double macrosigma = cross_sections::neutron_extrafirmsilicone_total_macro(nKE);
// dtnip = pdf::exponential(macrosigma); 
// step = gGeoManager->GetStep();

// if(step*inLV.P()>dtnip){
//   gGeoManager->SetStep(dtnip/inLV.P()); //only step to next interaction point
//   gGeoManager->Step();
//   if(global::myRnd.Rndm()<cross_sections::neutron_extrafirmsilicone_hydrogen_macro(nKE) / cross_sections::neutron_extrafirmsilicone_total_macro(nKE)){
//     sim::ElasticScatteringISO(inLV,M_P);
//   }
//   else{
//     sim::ElasticScatteringISO(inLV,M_C);
//   }
//   gGeoManager->SetCurrentDirection(inLV.Px(),inLV.Py(),inLV.Pz());
//   steptime += fabs(dtnip*M_N/(inLV.P()*300));
//   stepcounter++;
//   continue;
// }
// else{//step all the way to boundary and cross
//   gGeoManager->SetStep(step+(0.01/inLV.P())); //cross by 0.1 mm
//   gGeoManager->Step(kFALSE);
//   steptime += fabs(step*M_N/300); //dR = step *inLV.P()
//   stepcounter++;
//   continue;
// }
//     }



//     if(nodenum == kABSORBER){
// dtnip = GetAbsorberInteractionPoint(nKE);
// gGeoManager->FindNextBoundary();
// step = gGeoManager->GetStep();
// if(step*inLV.P()>dtnip){
//   gGeoManager->SetStep(dtnip/inLV.P()); //only step to next interaction point
//   gGeoManager->Step();
//   Absorber_hit(inLV);
//   gGeoManager->SetCurrentDirection(inLV.Px(),inLV.Py(),inLV.Pz());
//   steptime += fabs(dtnip*M_N/(inLV.P()*300));
//   stepcounter++;
//   continue;
// }
// else{//step all the way to boundary and cross
//   gGeoManager->SetStep(step+(0.01/inLV.P())); //cross by 0.1 mm
//   gGeoManager->Step(kFALSE);
//   steptime += fabs(step*M_N/300); //dR = step *inLV.P()
//   stepcounter++;
//   continue;
// }
//     }

//   }while(!gGeoManager->IsOutside());
//   return 0;
// }


// double AluminiumHousing::Al_hit(TLorentzVector& inLV){

//   //Create the Neut+Target LV
//   TLorentzVector neut_LVcopy(0,0,inLV.P(),inLV.E());

//   TLorentzVector Target(0.,0.,0.,M_AL);
//   TLorentzVector Before(Target + neut_LVcopy);
  
//   //Boost into CM frame
//   TVector3 boostv = (-1)*Before.BoostVector();
//   neut_LVcopy.Boost(boostv);
  
//   //increment Aluminium hit counter and see how much energy gets lost
//   if(fgNeutHit == 0){ //i'm only interested in # of 
//     //aluminium hits before I hit a neutron detector
//     fAluminiumCounter++;
//   }
//   double nKE = inLV.E()-inLV.M();
//   double theta(0),phi(0);
  
//   if(nKE<=0.35){
//     theta = acos(-1. + 2.*global::myRnd.Rndm()); //isotropic
//   }
//   else if(nKE>0.35&&nKE<=0.55){
//     theta = angdis[0].GetTheta();
//   }
//   else if(nKE>0.55&&nKE<=0.9){
//     theta = angdis[1].GetTheta();
//   }
//   else if(nKE>0.9&&nKE<=2.25){
//     theta = angdis[2].GetTheta();
//   }
//   else{
//     theta = angdis[3].GetTheta();
//   }
  
//   phi = 2.* 3.14 * global::myRnd.Rndm(); 

  
//   //Set CM LV and boost back to lab frame
//   neut_LVcopy.SetTheta(theta);
//   neut_LVcopy.SetPhi(phi);
//   neut_LVcopy.Boost(-1*boostv);//boost back to lab frame


//   //Rotate vector back from z-axis
//   neut_LVcopy.RotateY(inLV.Theta());
//   neut_LVcopy.RotateZ(inLV.Phi());
  

//   //Get Energy Loss
//   fELost += inLV.E()- neut_LVcopy.E();//add what neut loses
  
//   //return the new E/angle to the original neut Lorentz Vector
//   inLV = neut_LVcopy;

//   return theta;  
// }




// double AluminiumHousing::Absorber_hit(TLorentzVector& inLV){

//   //Create the Neut+Target LV
//   TLorentzVector neut_LVcopy(0,0,inLV.P(),inLV.E());
//   TLorentzVector Target(0.,0.,0.,M_CU);
//   TLorentzVector Before(Target + neut_LVcopy);
  
//   //Boost into CM frame
//   TVector3 boostv = (-1)*Before.BoostVector();
//   neut_LVcopy.Boost(boostv);
  
//   //increment Aluminium hit counter and see how much energy gets lost
//   fAbsorberCounter++;
  
//   //Calculate En CM + random CM_Angles
//   double theta = acos(-1. + 2.*global::myRnd.Rndm());
//   double phi = 2.* 3.14 * global::myRnd.Rndm(); //isotropic CM
  
//   //Set CM LV and boost back to lab frame
//   neut_LVcopy.SetTheta(theta);
//   neut_LVcopy.SetPhi(phi);
//   neut_LVcopy.Boost(-1*boostv);//boost back to lab frame
  
//   neut_LVcopy.SetTheta(neut_LVcopy.Theta()+inLV.Theta());
//   neut_LVcopy.SetPhi(neut_LVcopy.Phi()+inLV.Phi());
  
//   //Get Energy Loss
//   fELost += inLV.E()- neut_LVcopy.E();//add what neut loses
  
//   //return the new E/angle to the original neut Lorentz Vector
//   inLV = neut_LVcopy;

//   return theta;  
// }




// //find where on the AluminumHousing we hit.
// bool AluminiumHousing::inDet(const TVector3& v, const TVector3& beamspot){
//   //first see if it intersects the plane of the detector
//   TVector3 normv(0,0,fPos.Z()+kMidThick+kFrontThick);
//   double udotnormv = (v.Unit()).Dot(normv);
//   if(udotnormv <= 0.)
//     return false;
  
//   //Find distance from target to interesection point of detector plane 
//   //use line-plane intersection formula
//   TVector3 posvect(GetPosVect().X(),GetPosVect().Y(),fPos.Z() + kMidThick + kFrontThick);
//   double vdist = normv.Dot(posvect)/udotnormv;
//   if(vdist <= 0.)
//     return false;
  
//   //vector from target to interesection point of detector plane with
//   //magnitude equal to distance
//   TVector3 v_to_det(v);
//   v_to_det.SetMag(vdist);
  
//   //create vector from detector origin to interaction point on plane of detector
//   TVector3 ch_vect = v_to_det + beamspot - posvect;
  
//   //radial component from center of detector
//   double ch_vect_mag = ch_vect.Mag();
  
//   //see if it falls in the detector window
//   if((ch_vect_mag > fOuterrad)){
//     return false;
//   }
//   return true;
// }  




// int ReconstructTREL(NeutCollection&in,int&t_mult,double &t_first,int&det_first){
//   int n_tmult = 0;
//   double tfirst = 4096.0;
//   int detfirst = -1;
//   for(unsigned int i=0;i<in.size();i++){
//     if(in[i].T()>0){
// n_tmult++;
// if(in[i].T()<tfirst){
//   tfirst=in[i].T();
//   detfirst=i;
// }
//     }
//   }
//   if(tfirst<4096.0){    
//     //calculate TRel for all detectors(only important for coincidence data(ie source)
//     for(NeutCollectionRef it = in.begin(); it != in.end();it++){
// (*it).CalculateTRel(tfirst);
//     }
//     t_mult = n_tmult;
//     t_first = tfirst;
//     det_first = detfirst;
//     return 1;
  
//   }
//   t_mult = 0;
//   t_first = 0;
//   det_first = -1;
  
//   return 0;
// }



#endif
