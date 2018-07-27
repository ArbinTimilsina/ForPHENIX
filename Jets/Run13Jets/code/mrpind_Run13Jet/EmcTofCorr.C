#include <iostream>
#include <fstream>
#include <cmath>

#include <TOAD.h>
#include <TFile.h>
#include <TTree.h>
#include <EmcTofCorr.h>

using namespace std;

void ReadTOFMap(int fillnumber){

  //initialize the variable
  for(int i = 0; i < NARMSECTS_TOF; i++){
    for(int j = 0; j < NIYVAL; j++){
      for(int k = 0; k < NIZVAL; k++){
    
	tofmap[i][j][k] = 0;

      }//end k
    }//end j
  }//end i

  cout << "LOADING EMC TOF CORRECTION FILE" << endl;

  TOAD* tofmaploader = new TOAD("Run13Jet");
  string locationtofmap = tofmaploader->location("emcTOFmap/Run13pp510_EMC_TOF_Correction.root");
  
  TFile *tofmapf = new TFile(locationtofmap.c_str(), "read");
  
  if(!tofmapf->IsOpen())
    {
      cout << "Can not find " << locationtofmap.c_str() << "." << endl;
      delete tofmaploader;
      delete tofmapf;
      return;
    }

  TTree *Ttof = (TTree*)tofmapf->Get("T");
  if(Ttof==NULL){

    cout <<"Can not find T in tofmap root file." << endl;
    delete tofmaploader;
    delete tofmapf;

    return;
  }

  int fill;
  float tof[NARMSECTS_TOF][NIYVAL][NIZVAL];
  Ttof->SetBranchAddress("fillnumber", &fill);
  Ttof->SetBranchAddress("tof_correction", tof);
  
  int nentries = Ttof->GetEntries();
  bool chk=kFALSE;
  for(int i=0;i<nentries;i++){
    
    Ttof->GetEntry(i);
      
    if(fill==fillnumber){
  
      cout << "Correct TOF correction entry is found." << endl;
      cout << "Fill Number = " << fillnumber << endl;
      chk = kTRUE;
      break;
    }//endif(fillnumeber matching

  }//end i

  if(chk==kFALSE){

    cout << "Can not find correct TOF Correction entry." << endl;
    return;
  }
 
  for(int i=0;i<8;i++){
    for(int j=0;j<48;j++){
      for(int k=0;k<96;k++){
	
	tofmap[i][j][k] = tof[i][j][k];
      }
    }
  }

  delete tofmaploader;
  delete tofmapf;
  //  delete Ttof;
    
  return;

}//end ReadTOFMap

float TofCorrection(int armsect, int iy, int iz)
{
  
  float tofcorrection = 0;
  
  tofcorrection = tofmap[armsect][iy][iz];

  if(tofcorrection < -999) return 0;
  if(isinf(tofcorrection)) return 0;
  if(isnan(tofcorrection)) return 0;

  return tofcorrection;

}//end TofCorrection()
