#ifndef __IonChamber_cpp
#define __IonChamber_cpp

#include "IonChamber.hpp"

IonChamber::IonChamber(){};
IonChamber::~IonChamber(){};

void IonChamber::init(std::string name){
    _Name=name;
    HitPos.SetXYZ(0,0,0);
    _ZPos=0.;//mm
    _YPos=0.;
    _XPos=0.;
    _WireSeparation = 3.;//mm

    _RotVect.SetXYZ(0.,0.,0.);
    _NormVect.SetXYZ(0.,0.,0.);

    _IrisOpening=_WindowMax;

    //maximum is the sum of 3 grids with 4096 chs 
    xcluster.SetELimits(0.,4096.*3.); 
    ycluster.SetELimits(0.,4096.*3.);

    ESeg.Initialize(name+".ESeg",1);
    DESeg.Initialize(name+".DESeg",1);
    xgrid.Initialize(name+".xgrid",32);
    ygrid.Initialize(name+".ygrid",32);
    xcluster.Initialize(name+".xcluster",32);
    ycluster.Initialize(name+".ycluster",32);
}


void IonChamber::Reset(){
    HitPos.SetXYZ(0.,0.,0.);
    ESeg.Reset();
    DESeg.Reset();
    xgrid.Reset();
    ygrid.Reset();
    xcluster.Reset();
    ycluster.Reset();
}

void IonChamber::SetCalibrations(VariableMap& VarMap){
    ESeg.SetCalibrations(VarMap);
    DESeg.SetCalibrations(VarMap);
    xgrid.SetCalibrations(VarMap);
    ygrid.SetCalibrations(VarMap);
    xcluster.SetCalibrations(VarMap);
    ycluster.SetCalibrations(VarMap);

    VarMap.GetParam(Name()+".p2",_p2);
    VarMap.GetParam(Name()+".p1",_p1);
    VarMap.GetParam(Name()+".p0",_p0);

    VarMap.GetParam(Name()+".xpos",_XPos);
    VarMap.GetParam(Name()+".ypos",_YPos);
    VarMap.GetParam(Name()+".zpos",_ZPos);
    VarMap.GetParam(Name()+".wire_dist",_WireSeparation);

    double tempx=0,tempy=0,tempz=0;
    VarMap.GetParam(Name()+".xrot",tempx);
    VarMap.GetParam(Name()+".yrot",tempy);
    VarMap.GetParam(Name()+".zrot",tempz);

    _RotVect.SetXYZ(tempx,tempy,tempz);

    VarMap.GetParam(Name()+".irisopening",_IrisOpening);

    if(_IrisOpening > _WindowMax){
        printf("Iris opening set to larger than window maximum, using max value ");
        _IrisOpening = _WindowMax;
    }

    CalcNormVect();
}




void IonChamber::Print(){
    printf("-------------------------------------------------\n");
    printf("Detector information for %s \n",Name().c_str());
    printf("X position = %f \n",_XPos);
    printf("Y position = %f \n",_YPos);
    printf("-------------------------------------------------\n");
}




double IonChamber::SumE_X() const{
    double e = 0;
    for(int i=0;i<xgrid.Mult();i++){
        e += xgrid.E(i);
    }
    return e; 
}
double IonChamber::SumE_Y() const{
    double e = 0;
    for(int i=0;i<ygrid.Mult();i++){
        e += ygrid.E(i);
    }
    return e; 
}




//----------------------------------------------------------------------
// Pos_X only works if we are certain that there are not pile up events
// in the ion chamber since it's the weighted average of ALL grid wires.
//----------------------------------------------------------------------
double IonChamber::PosX(){          
    double pos = 0;
    double chA = 0;
    double chB = 0;
    for(int i=0;i<xgrid.Mult();i++){
        chA += xgrid.E(i)*xgrid.Ch(i);
        chB += xgrid.E(i);
    }
    if(chB>0){
        pos = chA/chB;
    }
    return pos;  
}
double IonChamber::PosY(){
    double pos = 0;
    double chA = 0;
    double chB = 0;
    for(int i=0;i<ygrid.Mult();i++){
        chA += ygrid.E(i)*ygrid.Ch(i);
        chB += ygrid.E(i);
    }
    if(chB>0){
        pos = chA/chB;
    }
    return pos;
}                                               



//--------------------------------------------------------
//For this method, xgrid and ygrid must sorted by energy
//--------------------------------------------------------
double IonChamber::PosXAdj(){
    bool check = false;
    double ch = 0;
    double weighted_e = 0;
    double total_e = 0;

    weighted_e = xgrid.E(0)*xgrid.Ch(0);
    total_e = xgrid.E(0);

    for(int i=1;i<xgrid.Mult();i++){
        if (((int)xgrid.Ch(i) == ((int)xgrid.Ch(0) + 1)) || ((int)xgrid.Ch(i) == ((int)xgrid.Ch(0) - 1))){
            weighted_e += xgrid.E(i)*xgrid.Ch(i);
            total_e += xgrid.E(i);      
            check = true;
        }
    }

    if(check && total_e > 0){
        ch = weighted_e/total_e;
    }else{
        ch = xgrid.Ch(0) + ZeroToOne->Rndm() - 0.5;
    }

    return ch; 
}  
double IonChamber::PosYAdj(){
    bool check = false;
    double ch = 0;
    double weighted_e = 0;
    double total_e = 0;

    weighted_e = ygrid.E(0)*ygrid.Ch(0);
    total_e = ygrid.E(0);

    for(int i=1;i<ygrid.Mult();i++){
        if( ((int)ygrid.Ch(i) == ((int)ygrid.Ch(0) + 1)) || ((int)ygrid.Ch(i) == ((int)ygrid.Ch(0) - 1)) ){
            weighted_e += ygrid.E(i)*ygrid.Ch(i);
            total_e += ygrid.E(i);      
            check = true;
        }
    }

    if(check && total_e > 0){
        ch = weighted_e/total_e;
    }else{
        ch = ygrid.Ch(0) + ZeroToOne->Rndm() - 0.5;
    }
    return ch; 
}                                               


//--------------------------------------------------------
// reconstruct hit position and fill HitPos vector
//--------------------------------------------------------
void IonChamber::ReconstructHitPos(double chx,double chy){
    // if(xgrid.Ch() >= 0 && ygrid.Ch() >= 0){ 
    if(chx >= 0 && chy >= 0){ 
        HitPos.SetX( (chx + ZeroToOne->Rndm() + _XPos)* _WireSeparation);
        HitPos.SetY( (chy + ZeroToOne->Rndm() + _YPos)* _WireSeparation);

        // HitPos.SetX( (chx + _XPos)* _WireSeparation);
        // HitPos.SetY( (chy + _YPos)* _WireSeparation);

        // HitPos.SetX( (xgrid.Ch() + ZeroToOne->Rndm() + _XPos)* _WireSeparation);
        // HitPos.SetY( (ygrid.Ch() + ZeroToOne->Rndm() + _YPos)* _WireSeparation);

        // HitPos.SetX((PosX() - 16.)* _WireSeparation);
        // HitPos.SetY((PosY() - 16.)* _WireSeparation);
        // HitPos.SetX((PosXAdj() - 16.) * _WireSeparation);
        // HitPos.SetY((PosYAdj() - 16.) * _WireSeparation);

        // if(xgrid.Mult() <= 1){
        //   HitPos.SetX( (xgrid.Ch(0) + (ZeroToOne->Rndm()-0.5)/2 + _XPos)* _WireSeparation);
        // }else if( xgrid.Ch(0)-xgrid.Ch(1) > 0){ //Event is on the left of wire
        //   HitPos.SetX((xgrid.Ch(0) - ZeroToOne->Rndm() + _XPos)* _WireSeparation);
        // }else if( xgrid.Ch(0)-xgrid.Ch(1) < 0){ //Event is on the right of wire
        //   HitPos.SetX( (xgrid.Ch(0) + ZeroToOne->Rndm() + _XPos)* _WireSeparation);
        // }
        // if(ygrid.Mult() <= 1){
        //   HitPos.SetY( (ygrid.Ch(0) + (ZeroToOne->Rndm()-0.5)/2 + _YPos)* _WireSeparation);
        // }else if( ygrid.Ch(0)-ygrid.Ch(1) > 0){ //Event is on the left of wire
        //   HitPos.SetY((ygrid.Ch(0) - ZeroToOne->Rndm() + _YPos)* _WireSeparation);
        // }else if( ygrid.Ch(0)-ygrid.Ch(1) < 0){ //Event is on the right of wire
        //   HitPos.SetY( (ygrid.Ch(0) + ZeroToOne->Rndm() + _YPos)* _WireSeparation);
        // }

        HitPos.SetZ(_ZPos);

        HitPos.RotateX(_RotVect[0]*TMath::DegToRad());
        HitPos.RotateY(_RotVect[1]*TMath::DegToRad());
        HitPos.RotateZ(_RotVect[2]*TMath::DegToRad());

    }else{
        HitPos.SetXYZ(-500,-500,-500);
        HitPos.SetMagThetaPhi(0,0,0);
    }

}






void IonChamber::CalcNormVect(){
    TVector3 resultv = TVector3(_XPos,_YPos,_ZPos).Unit();
    if(_RotVect[0]){
        resultv.RotateX(_RotVect[0]*TMath::DegToRad());
    }
    if(_RotVect[1]){
        resultv.RotateY(_RotVect[1]*TMath::DegToRad());
    }
    if(_RotVect[2]){
        resultv.RotateZ(_RotVect[2]*TMath::DegToRad());
    }
    _NormVect = resultv;
}

bool IonChamber::InDet(const TVector3& v, const TVector3& beamspot){
    //first see if it intersects the plane of the detector
    double udotnormv = (v.Unit()).Dot(_NormVect);
    if(udotnormv <= 0.)
        return false;

    //Find distance from target to interesection point of detector plane 
    //use line-plane intersection formula
    TVector3 posvect = GetPosVect();
    double vdist = _NormVect.Dot(posvect)/udotnormv;
    if(vdist <= 0.)
        return false;

    //vector from target to interesection point of detector plane with
    //magnitude equal to distance
    TVector3 v_to_det(v);
    v_to_det.SetMag(vdist);

    //create vector from detector origin to interaction point on plane of detector
    TVector3 ch_vect = v_to_det + beamspot - posvect;

    //radial component from center of detector
    double ch_vect_mag = ch_vect.Mag();

    //see if it falls in the detector window
    if((ch_vect_mag > _IrisOpening)){
        return false;
    }
    return true;
}



bool IonChamber::Vect_to_ch(const TVector3& v,double& cx,double& cy, const TVector3& beamspot){
    //first see if it intersects the plane of the detector
    double udotnormv = (v.Unit()).Dot(_NormVect);
    if(udotnormv <= 0.)
        return false;

    //Find distance from target to interesection point of detector plane 
    //use line-plane intersection formula
    TVector3 posvect = GetPosVect();
    double vdist = _NormVect.Dot(posvect)/udotnormv;
    if(vdist <= 0.)
        return false;

    //vector from target to interesection point of detector plane with
    //magnitude equal to distance
    TVector3 v_to_det(v);
    v_to_det.SetMag(vdist);

    //create vector from detector origin to interaction point on plane of detector
    TVector3 ch_vect = v_to_det + beamspot- posvect;

    //radial component from center of detector
    double ch_vect_mag = ch_vect.Mag();

    //see if it falls in the detector window
    if((ch_vect_mag > _IrisOpening)){
        return false;
    }

    //find the wires
    //TMath::Nint rounds to nearest integer
    cx = TMath::Nint((ch_vect.X() / _WireSeparation) + (double)xgrid.NumOfCh()/2);
    cy = TMath::Nint((ch_vect.Y() / _WireSeparation) + (double)ygrid.NumOfCh()/2);

    if(cx<0 || cx >=  xgrid.NumOfCh()){
        return false;
    }
    if(cy<0 || cy>= ygrid.NumOfCh()){
        return false;
    }
    return true;
}




/*************************************************************************************
  Reconstruction of clusters 
 ************************************************************************************/
void IonChamber::ReconstructCluster(Detector& grid,Detector& cluster)
{

  // if there is no multiplicity add the single channel into the cluster and exit the Reconstruction function
  if(grid.Mult() < 2 && grid.Mult() > 0)
  {
    cluster.InsertHit(grid.E(),0,grid.Ch());
    cluster.SortByEnergy();
    return;
  } else { // find neighbors
    for(int i=0; i < grid.Mult()-1; i++) // loop through 0 -> end-1
    {
      for(int j=i+1; j < grid.Mult(); j++) // loop through the i+1 -> end 
      {
        if( i != j && (TMath::Abs(grid.Ch(i)-grid.Ch(j)) < 2  && TMath::Abs(grid.Ch(i)-grid.Ch(j)) > 0) )
        {
          cluster.InsertHit(grid.E(i) + grid.E(j), 0, (grid.Ch(i)*grid.E(i) + grid.Ch(j)*grid.E(j))/(grid.E(i) + grid.E(j)));
        }
      }
    }
  }

  cluster.SortByEnergy();
}

/* void IonChamber::ReconstructCluster(Detector& grid,Detector& cluster) */
/* { */
/*     // if there is no multiplicity add the single channel into the cluster and exit the Reconstruction function */
/*     if(grid.Mult() < 2 && grid.Mult() > 0) */
/*     { */
/*         cluster.InsertHit(grid.E(),0,grid.Ch()); */
/*         cluster.SortByEnergy(); */
/*         return; */
/*     } */

/*     // get the size of the Multiplicity for the size of the raw arrays and fill them with zero */
/*     std::vector<double> e(grid.Mult(),0.0); */ 
/*     std::vector<double> ch(grid.Mult(),0.0); */

/*     // fill the channel and energy vector */
/*     for(int i=0;i<grid.Mult();i++) */
/*     { */
/*         e[i]=grid.E(i); */
/*         ch[i]=grid.Ch(i); */
/*     } */

/*     // set the pointers of the two iterators we will use to step through */
/*     std::vector<double>::iterator e1_it = e.begin(); */
/*     std::vector<double>::iterator ch1_it = ch.begin(); */
/*     std::vector<double>::iterator e2_it = e.begin() + 1 ; // point to the second member of ecluster */
/*     std::vector<double>::iterator ch2_it = ch.begin() + 1; // point to the second member of chcluster */

/*     while(e2_it < e.end()) // this is the condition to look through the array until reaching the end */
/*     { */
/*         while(e2_it != e1_it) // this condition enforces that the iterators are never the same! */
/*         { */
/*             // is there a neighboring channel? */
/*             if(TMath::Abs(*ch2_it - *ch1_it) < 2 && TMath::Abs(*ch2_it - *ch1_it) > 0) // look for channels that are only 1 away */
/*             { */
/*                 *ch1_it = ((*ch1_it)*(*e1_it) + (*ch2_it)*(*e2_it)) / (*e1_it + *e2_it); // calc an energy weighted location */ 
/*                 *e1_it += *e2_it; // sum the total energies */
/*                 e.erase(e2_it); */
/*                 ch.erase(ch2_it); */
/*                 e2_it--; */
/*                 ch2_it--; */
/*                 break; // break out of the while(e2_it != e1_it) */
/*             } */
/*             else // if no neighboring channel move to the next index */
/*             { */
/*                 e1_it++; */
/*                 ch1_it++; */
/*             } */
/*         } */
/*         e1_it = e.begin(); // put the e1 back into the first location */ 
/*         ch1_it = ch.begin(); */
/*         e2_it++; // set to the next index */
/*         ch2_it++; */
/*     } */

/*     //insert hits into xcluster detector */
/*     ch1_it = ch.begin(); */
/*     for(auto it = e.begin(); it < e.end(); ++it){ */
/*         // for(std::vector<double>::iterator it = e.begin(); it < e.end(); ++it){ */
/*         cluster.InsertHit(*it,0,*ch1_it); */
/*         ++ch1_it; */
/*     } */

/*     cluster.SortByEnergy(); */
/*     } */




    // void IonChamber::ReconstructXClusters(){
    //   if(xgrid.Mult() <= 0){
    //     return;
    //   }

    //   std::vector<double> ecluster(xgrid.Mult(),0.0); //get the size of the raw arrays
    //   std::vector<double> chcluster(xgrid.Mult(),0.0);

    //   for(int i=0;i<xgrid.Mult();i++){
    //     ecluster[i] = xgrid.E(i);
    //     chcluster[i] = xgrid.Ch(i);
    //   }

    //   std::vector<double>::iterator e2_it = ecluster.begin() + 1 ; //point to the second member of ecluster
    //   std::vector<double>::iterator ch2_it = chcluster.begin() + 1;//point to the second member of chcluster
    //   std::vector<double>::iterator ch1_it = chcluster.begin();
    //   std::vector<double>::iterator e1_it = ecluster.begin();

    //   while(e2_it < ecluster.end()){
    //     while(e2_it != e1_it){
    //       if(TMath::Abs(*ch2_it - *ch1_it) < 2){
    //         *ch1_it = ((*ch1_it)*(*e1_it) + (*ch2_it)*(*e2_it)) / (*e1_it + *e2_it);
    //         *e1_it += *e2_it;
    //         ecluster.erase(e2_it);
    //         chcluster.erase(ch2_it);
    //         e2_it--;
    //         ch2_it--;
    //         break;
    //       }
    //       else{
    //         e1_it++;
    //         ch1_it++;
    //       }
    //     }
    //     e1_it = ecluster.begin();
    //     ch1_it = chcluster.begin();
    //     e2_it++;
    //     ch2_it++;
    //   }

    //   //insert hits into xcluster detector
    //   ch1_it = chcluster.begin();
    //   for(std::vector<double>::iterator it = ecluster.begin(); it<ecluster.end(); ++it){
    //     xcluster.InsertHit(*it,0,*ch1_it);
    //     ++ch1_it;
    //   }
    //   xcluster.SortByEnergy();
    // }

    // void IonChamber::ReconstructYClusters(){
    //   if(ygrid.Mult()<=0){
    //     return;
    //   }
    //   std::vector<double> ecluster((int)ygrid.Mult(),0); //get the size of the raw arrays
    //   std::vector<double> chcluster((int)ygrid.Mult(),0);

    //   for(int i=0;i<ygrid.Mult();i++){  
    //     ecluster[i] = ygrid.E(i);
    //     chcluster[i] = ygrid.Ch(i);
    //   }

    //   std::vector<double>::iterator e2_it = ecluster.begin() + 1; ;//point to the second member of ecluster
    //   std::vector<double>::iterator ch2_it = chcluster.begin() + 1;//point to the second member of chcluster
    //   std::vector<double>::iterator ch1_it = chcluster.begin();//point to first element
    //   std::vector<double>::iterator e1_it = ecluster.begin();//point to first element

    //   while(e2_it < ecluster.end()){
    //     while(e2_it != e1_it){ //look at all previous entries
    //       if(TMath::Abs(*ch2_it - *ch1_it) < 2){ //compare the channel of all previous entries to this channel
    //       	*ch1_it = ((*ch1_it)*(*e1_it) + (*ch2_it)*(*e2_it)) / (*e1_it + *e2_it);
    //       	*e1_it += *e2_it;
    //       	ecluster.erase(e2_it);
    //       	chcluster.erase(ch2_it);
    //       	e2_it--;
    //       	ch2_it--;
    // 	      break;
    //       }
    //       else{
    //       	e1_it++;
    //       	ch1_it++;
    //       }
    //     }
    //     e1_it = ecluster.begin();
    //     ch1_it = chcluster.begin();
    //     e2_it++;
    //     ch2_it++;
    //   }

    //   //insert hits into xcluster detector
    //   ch1_it = chcluster.begin();
    //   for(std::vector<double>::iterator it = ecluster.begin();it<ecluster.end();++it){
    //     ycluster.InsertHit(*it,0,*ch1_it);
    //     ch1_it++;
    //   }

    //   ycluster.SortByEnergy();
    // }





#endif
