/***************************************************************
Class: SiliconDetector, SiliconCluster, SiliconArray

front silicon in telescope is si_a -> S2Detector
followed by the second, si_b -> S1Detector
****************************************************************/

#ifndef __SiliconDetector__cpp
#define __SiliconDetector__cpp

#include <cmath>
#include "SiliconDetector.hpp"

void SiliconDetector::init(std::string n, const int& fnum,const int& bnum){
  name=n;
  frontNum=fnum;
  backNum=bnum;

  Front.Initialize(name+".Front",frontNum);
  Back.Initialize(name+".Back",backNum);
  FrontCluster.Initialize(name+".FrontCluster",frontNum);
  BackCluster.Initialize(name+".BackCluster",backNum);
  
  Ch_front.resize(frontNum);
  Ch_back.resize(backNum);
  EMatch.resize(frontNum);
  TMatch.resize(frontNum);
  Pos.resize(frontNum);

  // outerRad = S2OUTERRAD; // J. Baker only had S2 detectors!
  // innerRad = S2INNERRAD; // Remeber to changes these back for Ne-19
  // if(name == "S2"){
  //   outerRad = S2OUTERRAD;
  //   innerRad = S2INNERRAD;
  // }else if(name == "S1"){
  //   outerRad = S1OUTERRAD;
  //   innerRad = S1INNERRAD;
  // }else{
  //   std::cerr << "Cannot set silicon parameters for:" << name << std::endl;
  //   exit(EXIT_FAILURE);
  // }  
  // ringPitch = abs(outerRad - innerRad)/16.0;
  // deltaPhi = 2.0*TMath::Pi()/16.0;

  normVect.SetXYZ(0.,0.,0.);
  posVect.SetXYZ(0.,0.,0.);
}


void SiliconDetector::Reset(){
  Front.Reset(); Back.Reset();
  FrontCluster.Reset(); BackCluster.Reset();
  Ch_front.clear(); Ch_back.clear();
  EMatch.clear(); TMatch.clear();
  Pos.clear();
}

void SiliconDetector::Print(){
  printf("-------------------------------------------------\n");
  printf("Detector information for %s \n",Name().c_str());
  printf("Match Energy From Back Ch = %s \n", UseBackEnergy ? "true" : "false");
  printf("Match Epsilon = %f \n",matchEpsilon);
  printf("Match Delta = %f \n",matchDelta);
  printf("Cluster Front = %s \n",f_cluster ? "true" : "false");
  printf("Cluster Back = %s \n",b_cluster ? "true" : "false");
  printf("Inner radius = %f \n",innerRad);
  printf("Outer radius = %f \n",outerRad);
  printf("Ring Pitch = %f \n",ringPitch);
  printf("Delta Phi = %f \n",deltaPhi);
  printf("Front Ch = %f \n",Front.Ch());
  printf("Back Ch = %f \n",Back.Ch());
  printf("Theta = %f \n",Theta()*TMath::RadToDeg());
  printf("Phi = %f \n",Phi()*TMath::RadToDeg());
  printf("-------------------------------------------------\n");
}


void SiliconDetector::SetCalibrations(VariableMap& VarMap){  
  Front.SetCalibrations(VarMap);
  Back.SetCalibrations(VarMap);
  // FrontCluster.SetCalibrations(VarMap); //Not being used right now
  // BackCluster.SetCalibrations(VarMap);  //Not being used right now

  VarMap.GetParam(Name()+".Rin",innerRad);
  VarMap.GetParam(Name()+".Rout",outerRad);
  ringPitch = abs(outerRad - innerRad)/16.0;
  deltaPhi = 2.0*TMath::Pi()/16.0;

  //Calibratrions for cluster reconstruction
  VarMap.GetParam(Name()+".UseBackEnergy",UseBackEnergy);
  VarMap.GetParam(Name()+".MatchEpsilon",matchEpsilon);
  VarMap.GetParam(Name()+".MatchDelta",matchDelta);
  VarMap.GetParam(Name()+".Cluster_Front",f_cluster);
  VarMap.GetParam(Name()+".Cluster_Back",b_cluster);
  
  double tempx=0.,tempy=0.,tempz=0.;
  VarMap.GetParam(Name()+".xpos",tempx);
  VarMap.GetParam(Name()+".ypos",tempy);
  VarMap.GetParam(Name()+".zpos",tempz);
  posVect.SetXYZ(tempx,tempy,tempz);
  
  VarMap.GetParam(Name()+".xrot",rot[0]);
  VarMap.GetParam(Name()+".yrot",rot[1]);
  VarMap.GetParam(Name()+".zrot",rot[2]);
  
  CalcNormVect();
}



void SiliconDetector::CalcNormVect(){
  // TVector3 resultv = (_shiftVect + _posVect).Unit();
  TVector3 resultv =  posVect.Unit();
  resultv.RotateX(rot[0]*TMath::DegToRad());
  resultv.RotateY(rot[1]*TMath::DegToRad());
  resultv.RotateZ(rot[2]*TMath::DegToRad());  
  normVect = resultv;
}


double SiliconDetector::ThetaMin(){
  return TMath::ATan(innerRad/posVect.Z());
}
double SiliconDetector::ThetaMax(){
  return TMath::ATan(outerRad/posVect.Z());
}

/* This function if for simulation and efficiency calculation */
bool SiliconDetector::InDet(const TVector3& sim_vect, const TVector3& beamspot){
  //first see if it intersects the plane of the detector
  double udot_normvect = (sim_vect.Unit()).Dot(normVect);
  if(udot_normvect <= 0.){
    return false;
  }
  //Find distance from target to interesection point of detector plane 
  //use line-plane intersection formula
  double vdist = normVect.Dot(posVect)/udot_normvect;
  if(vdist <= 0.){
    return false;
  }
  //vector from target to interesection point of detector plane with
  //magnitude equal to distance
  TVector3 v_to_det(sim_vect);
  v_to_det.SetMag(vdist);
  //create vector from detector origin to interaction point on plane of detector
  TVector3 ch_vect = v_to_det + beamspot - posVect;
  //radial component from center of detector
  double ch_vect_mag = ch_vect.Mag();
  //see if it falls in the detector window
  if((ch_vect_mag > outerRad) || (ch_vect_mag < innerRad)){
    return false;
  }
  return true;
}




/*****************************************************************************************
This Class Ch2Vect converts the cf (front channel) and cb (back channel) to TVectors named
resultVect using a Rndm number between channels
******************************************************************************************/
TVector3 SiliconDetector::Ch2Vect(double ch_front, double ch_back){
  double s = 0.0, phi = 0.0;
  TVector3 resultVect(0.,0.,0.);  

  s = innerRad + (ch_front + Rndm->Uniform(-0.5,0.5))*ringPitch; // et commented for jess' data
  phi = (ch_back + Rndm->Uniform(0,1.))*deltaPhi;
  // s = innerRad + ch_front*ringPitch;
  // phi = ch_back*deltaPhi;

  resultVect.SetMagThetaPhi(s, TMath::Pi()/2.0, phi); //Everything in RADIANS!
  resultVect.RotateX(rot[0]*TMath::DegToRad());
  resultVect.RotateY(rot[1]*TMath::DegToRad());
  resultVect.RotateZ(rot[2]*TMath::DegToRad());

  // Pos[0] = resultVect + posVect;
  return resultVect + posVect;
}



/******************************************************************
The ReconstuctClusters class take a detector S1 or S2 and the 
channels associated with such detector and combines them into a 
cluster (grouping)

The S2/S1Detector should already be sorted by energy before calling
the reconstruction clusters
******************************************************************/
void SiliconDetector::ReconstructClusters(SiliconDetector& Det){
  int i=0,j=0;
  double e_cluster = 0.0, t_cluster = 0.0, ch_cluster = 0.0, e_tot = 0.0;
  // bool front_status = false, back_status = false;

  if(Det.Back.Mult() < 2 && Det.Back.Mult() > 0)
  {
    BackCluster.InsertHit(Det.Back.E(),Det.Back.T(),Det.Back.Ch());
  } else { // find neighbors
    for(int i=0; i < Det.Back.Mult()-1; i++) // loop through 0 -> end-1
    {
      for(int j=i+1; j < Det.Back.Mult(); j++) // loop through the i+1 -> end 
      {
        if( i != j && (TMath::Abs(Det.Back.Ch(i)-Det.Back.Ch(j)) < 2  && TMath::Abs(Det.Back.Ch(i)-Det.Back.Ch(j)) > 0) )
        {
          BackCluster.InsertHit(Det.Back.E(i) + Det.Back.E(j), (Det.Back.T(j+1) + Det.Back.T(j))/2.0, (Det.Back.Ch(i)*Det.Back.E(i) + Det.Back.Ch(j)*Det.Back.E(j))/(Det.Back.E(i) + Det.Back.E(j)));
        }
      }
    }
  }
  BackCluster.SortByEnergy();

  if(Det.Front.Mult() < 2 && Det.Front.Mult() > 0)
  {
    FrontCluster.InsertHit(Det.Front.E(),Det.Front.T(),Det.Front.Ch());
    // FrontCluster.SortByEnergy();
    // return;
  } else { // find neighbors
    for(int i=0; i < Det.Front.Mult()-1; i++) // loop through 0 -> end-1
    {
      for(int j=i+1; j < Det.Front.Mult(); j++) // loop through the i+1 -> end 
      {
        if( i != j && (TMath::Abs(Det.Front.Ch(i)-Det.Front.Ch(j)) < 2  && TMath::Abs(Det.Front.Ch(i)-Det.Front.Ch(j)) > 0) )
        {
          FrontCluster.InsertHit(Det.Front.E(i) + Det.Front.E(j), (Det.Front.T(j+1) + Det.Front.T(j))/2.0, (Det.Front.Ch(i)*Det.Front.E(i) + Det.Front.Ch(j)*Det.Front.E(j))/(Det.Front.E(i) + Det.Front.E(j)));
        }
      }
    }
  }
  FrontCluster.SortByEnergy();

  j=0; i=0; 

  double energy_sep = 0.0;
  double delta = 0.0;
  double back_e = 0.0;

  /* Front and Backclusters are already sorted by energy, 
  so i and j starts at zero i.e. the highest energy event. */
  while( (i < FrontCluster.Mult()) &&  (j < BackCluster.Mult()) ){
    energy_sep = TMath::Abs(FrontCluster.ERaw(i) - BackCluster.ERaw(j)); // Front - Back
    if( UseBackEnergy ){ 
      delta = matchEpsilon + matchDelta*BackCluster.ERaw(j);
    }else{
      delta = matchEpsilon + matchDelta*FrontCluster.ERaw(j);
    }
    // delta = matchEpsilon + matchDelta*BackCluster.ERaw(j);
    back_e = BackCluster.ERaw(j);

    if( TMath::Abs(energy_sep) <= delta ){ // we have a match
      double front_e = FrontCluster.ERaw(i);
      double front_ch = FrontCluster.ChRaw(i);
      double back_ch = BackCluster.ChRaw(j);
      double back_t = BackCluster.TRaw(j);

      Pos[j] = Ch2Vect(front_ch, back_ch);
      
      Ch_front[j] = front_ch;
      Ch_back[j] = back_ch;
      
      if( UseBackEnergy ){ 
        EMatch[j] = back_e;
        TMatch[j] = back_t;
      }else{
        EMatch[j] = front_e;
        TMatch[j] = back_t;
      }

      j++;

      // // sum up matched energies front and back
      // energyMatchFront += front_e;
      // energyMatchBack += back_e;

    //end found match
    }else if(energy_sep > delta){
      // Front energy is larger than back energy+delta ,so we won't ever match 
      // this one in the future, drop it, pick next front, keep the back
      // front_status = false; // flags "unmatched" 
      i++;
      continue;
    }else if(energy_sep < (-delta)){
      // back energy is larger than front energy-delta , 
      // so we won't ever match 
      // this one in the future, drop it, pick next back, keep the front
      // back_status = false; // flags "unmatched" 
      j++;
      continue;
    }

    i++;
    j++;
  }// end of while loop
  
  return;
}

// void SiliconDetector::ReconstructClusters(SiliconDetector& Det){
//   int i=0,j=0;
//   double e_cluster = 0.0, t_cluster = 0.0, ch_cluster = 0.0, e_tot = 0.0;
//   // bool front_status = false, back_status = false;

//   //------------------------------------------
//   // Look at the back wedges (using j index)!!!
//   //-----------------------------------------
//   for(j=0; j<Det.Back.Mult(); j++){
//     e_tot = e_tot + Det.Back.E(j); 

//     if(b_cluster && (i+1 < Det.Back.Mult()) && ( (Det.Back.Ch(i+1)==Det.Back.Ch(i)+1) || (Det.Back.Ch(i+1)==Det.Back.Ch(i)-1)))
//     {
//       e_cluster = Det.Back.E(j+1) + Det.Back.E(j);
//       t_cluster = (Det.Back.T(j+1) + Det.Back.T(j))/2.0;
//       // ch_cluster = (Det.Back.Ch(j+1) + Det.Back.Ch(j))/2.0;
//       ch_cluster = (Det.Back.Ch(j+1)*Det.Back.E(j+1) + Det.Back.Ch(j)*Det.Back.E(j))/e_cluster;

//       BackCluster.InsertHit(e_cluster,t_cluster,ch_cluster);
//       //we analyzed the j+1 hit along with the j, so 
//       j++;
//     }else{
//       BackCluster.InsertHit(Det.Back.E(i),Det.Back.T(i),Det.Back.Ch(i));
//     }
//   }


//   //-----------------------------------------------
//   // Now look at the front rings (using i index)!!!
//   //-----------------------------------------------
//   i=0;j=0;
//   e_cluster = 0.0; t_cluster = 0.0; ch_cluster = 0.0; e_tot = 0.0;

//   for(i=0; i<Det.Front.Mult(); i++){
//     e_tot = e_tot + Det.Front.E(i); 

//     if(f_cluster && (i+1<Det.Front.Mult()) && ( (Det.Front.Ch(i+1)==Det.Front.Ch(i)+1) || (Det.Front.Ch(i+1)==Det.Front.Ch(i)-1))){
//       e_cluster = Det.Front.E(i+1) + Det.Front.E(i);
//       t_cluster = (Det.Front.T(i+1) + Det.Front.T(i))/2.0;
//       // ch_cluster = (Det.Front.Ch(i+1) + Det.Front.Ch(i))/2.0;
//       ch_cluster = (Det.Front.Ch(j+1)*Det.Front.E(j+1) + Det.Front.Ch(j)*Det.Front.E(j))/e_cluster;

//       FrontCluster.InsertHit(e_cluster,t_cluster,ch_cluster);
//       //we analyzed the i+1 hit along with the i, so 
//       i++;
//     }else{
//       FrontCluster.InsertHit(Det.Front.E(i),Det.Front.T(i),Det.Front.Ch(i));
//     }
//   }
  
//   //---------------------------------------------------------------------------------
//   FrontCluster.SortByEnergy(); 
//   BackCluster.SortByEnergy(); 

//   // Pos[0] = Ch2Vect(FrontCluster.ChRaw(), BackCluster.ChRaw());
//   //---------------------------------------------------------------------------------  

//   j=0; i=0; 
//   // energyMatchFront=0.0;
//   // energyMatchBack=0.0;
  
//   double energy_sep = 0.0;
//   double delta = 0.0;
//   double back_e = 0.0;

//   /* Front and Backclusters are already sorted by energy, 
//   so i and j starts at zero i.e. the highest energy event. */
//   while( (i < FrontCluster.Mult()) &&  (j < BackCluster.Mult()) ){
//     energy_sep = FrontCluster.ERaw(i) - BackCluster.ERaw(j); // Front - Back
//     if( UseBackEnergy ){ 
//       delta = matchEpsilon + matchDelta*BackCluster.ERaw(j);
//     }else{
//       delta = matchEpsilon + matchDelta*FrontCluster.ERaw(j);
//     }
//     // delta = matchEpsilon + matchDelta*BackCluster.ERaw(j);
//     back_e = BackCluster.ERaw(j);

//     if( TMath::Abs(energy_sep) <= delta ){ // we have a match
//       double front_e = FrontCluster.ERaw(i);
//       double front_ch = FrontCluster.ChRaw(i);
//       double back_ch = BackCluster.ChRaw(j);
//       double back_t = BackCluster.TRaw(j);

//       Pos[j] = Ch2Vect(front_ch, back_ch);
      
//       Ch_front[j] = front_ch;
//       Ch_back[j] = back_ch;
      
//       if( UseBackEnergy ){ 
//         EMatch[j] = back_e;
//         TMatch[j] = back_t;
//       }else{
//         EMatch[j] = front_e;
//         TMatch[j] = back_t;
//       }

//       j++;

//       // // sum up matched energies front and back
//       // energyMatchFront += front_e;
//       // energyMatchBack += back_e;

//     //end found match
//     }else if(energy_sep > delta){
//       // Front energy is larger than back energy+delta ,so we won't ever match 
//       // this one in the future, drop it, pick next front, keep the back
//       // front_status = false; // flags "unmatched" 
//       i++;
//       continue;
//     }else if(energy_sep < (-delta)){
//       // back energy is larger than front energy-delta , 
//       // so we won't ever match 
//       // this one in the future, drop it, pick next back, keep the front
//       // back_status = false; // flags "unmatched" 
//       j++;
//       continue;
//     }

//     i++;
//     j++;
//   }// end of while loop
  

//   // frontMatchStat = front_status;
//   // backMatchStat = back_status;

//   return;
// }








/***************************************************************************
Silicon Array class
****************************************************************************/
void SiliconArray::init(std::string n){
  name=n;
  Si.Initialize(name,0);
}

// void SiliconArray::SetCalibrations(VariableMap& ){
//   VarMap->GetParam(Form("%s.thickness",Name()),fThickness); 
// }

// void SiliconArray::Reset(){
//   for(int i=0;i<fNumOfSi;i++){
//     fE_[i] = 0;
//     fT_[i] = 0;
//     fPos_[i].SetXYZ(0,0,0);
//   }
// }

// Double32_t SiliconArray::ERecoAB(){
//   double de = ERecoA();
//   if(de!=0){
//     return (de + E_B());
//   }
//   else{
//     return 0;
//   }
// }

// Double32_t SiliconArray::ERecoA(){
//   if(E_B()!= 0){
//     double path = 1/TMath::Cos(Theta_B());
//     double e = fELossTableSiA.GetLossLinear(path,E_B(),fThickness);
//     return e;
//   }  
//   else{
//     return 0;
//   }
// }




#endif
