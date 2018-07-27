#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "TAcceptParticle.h"
#include "TMath.h"

TAcceptParticle::TAcceptParticle()
{
  // MC info form dcgeom96.f
  RadiusDC = 224.0;
  PhiBotW = -0.589;
  PhiTopW = 0.982;
  PhiBotE = 3.731;
  PhiTopE = 2.160;
  
  Sl1 = 0.8;
  Sl2 = 0.1;
  
  fullTwoPi = false; 

}

int TAcceptParticle::acceptParticle(TLorentzVector *part, int charge, double vertex, 
				    bool &phiAccept, bool &zedAccept, bool &bigzedAccept, bool &deadArea )
{

  int val=0;
  double  zed, alpha, phi1, phi2; 

  double theta = part->Theta();
  double phi = part->Phi();
  double pt = part->P()*sin(theta); 
  double pi2 = TMath::Pi()/2.0;
  double pi = TMath::Pi();

  phiAccept = false; 
  zedAccept = false; 
  bigzedAccept = false; 
  deadArea = false; 

  // Simplified acceptance filter (adapted from S. Butsyk's fast acceptance simulator) 

  // zed at DCH assuming straight track

  zed  = vertex - RadiusDC*tan(theta-pi2);

  // Set phi to wrap at 3Pi/2

  if( phi<-pi2 ) phi+=2*pi;  
  if( phi>=3.0*pi2) phi-=2*pi; 

  // Dch zed geometry acceptance
  // Small dead area on side crossing, cut is approx abs(eta)<0.38 for vertex = 0
  if(!fullTwoPi){
    if ( fabs(zed)<0.5 ) deadArea = true; 
    if ( (fabs(zed)<90.0) && (fabs(zed)>0.5) ) zedAccept = true; 
  }
  else{
    if ( (fabs(zed)<90.0) ) zedAccept = true; 
  }

  // bigzedAccept is approximately abs(eta)<1.0 for vertex=0
  // (No dead area)
  if ( fabs(zed)<264.0 ) bigzedAccept = true; 
    
  if(pt>0.2){

      // Alpha vs. pt, zed parametrization from data
      alpha = -charge*(0.1015-((1.452e-6)*zed)+((2.356e-6)*zed*zed))/pt;      
      phi1 = phi-(alpha)*Sl1;
      phi2 = phi+(alpha)*Sl2;

      if(fullTwoPi){
    
        phiAccept = true; 

        if ( (phi>=pi2) && (phi<(3.0*pi/2.0)) ){  // Falls into East phi acceptance
          if ( zedAccept ) val = 1; 
        }
        else if ( (phi>=-pi2) && (phi<pi2) ){  // Falls into West phi acceptance
          if ( zedAccept ) val = 2; 
        }
        //else{
        //  cout << "WTF? phi = " << phi << endl; 
        //}

      }
      else{

        if ( (phi1>PhiTopE) && (phi1<PhiBotE) && (phi2>PhiTopE) && (phi2<PhiBotE) ){  // Falls into East phi acceptance
          phiAccept = true; 
          if ( zedAccept ) val = 1; 
        }
        else if ( (phi1>PhiBotW) && (phi1<PhiTopW) && (phi2>PhiBotW) && (phi2<PhiTopW) ){  // Falls into West phi acceptance
          phiAccept = true; 
          if ( zedAccept ) val = 2; 
        }

    }

  }
 
  return val;
}


