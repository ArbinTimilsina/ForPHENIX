//#include "CommonVar.h"
#include "EmcMatching.h"
#include <iostream>
#include <string>
#include <fstream>

#include <TOAD.h>
#include <TString.h>

using namespace std;

const double  cuthigh[NPT] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 10.0, 25.0};
const double  cutlow[NPT] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 10.0};

int GetPtBin(float pT)
{
  int kbin = -1;
  for (int m1 = 0; m1 < NPT; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

float GetPtBinCenter(int pt_bin)
{
  float pt_bin_center = (cutlow[pt_bin]+cuthigh[pt_bin])/2.;
  return pt_bin_center;
}

void LoadAllEMCFiles(){
  cout << "LOADING ALL THE FUCKING EMC FILES" << endl;
  //for sdphi and sdz
  LoadPoints();

  //for offset
  LoadMean();
  LoadSigma();

}//end LoadAllFiles()

void LoadPoints(){

  
  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++)
	{
	  for (int ized = 0; ized < NZED; ized++)
	    {
	      for (int icent = 0; icent < NCENT; icent++)
		{
		  
		  for (int ipt = 0; ipt < NPT; ipt++)
		    {
		      for(int i = 0; i < 3; i++)
			{
			  for(int j = 0; j < 3; j++)
			    {

			      mean_points_emc[icharge][iarmsect][ized][icent][ipt][i][j] = 0;
			      sigma_points_emc[icharge][iarmsect][ized][icent][ipt][i][j] = 0;

			    }//end j
			}//end i
		    }//end ipt
		}//end icent
	    }//end ized
	}//end iarmsect
    }//end icharge

  cout << "Attempting to execute LoadPoints " << endl;
  std::string dphi_dz_str[2] = {"dphi","dz"};
  std::string dPhi_dZ_str[2] = {"dPhi","dZ"};
  std::string Stage_str[3] = {"Initial","Intermediate","Final"};

  TOAD *toad_loader_mean = new TOAD("Run13Jet");
  TOAD *toad_loader_sigma = new TOAD("Run13Jet");

  for(int istage = 0; istage < 2; istage++){
    for(int idphi_dz = 0; idphi_dz < 2; idphi_dz++){

      if(istage == INTERMEDIATE || istage == FINAL){
	dPhi_dZ_str[DPHI] = "sdPhi";
	dPhi_dZ_str[DZ] = "sdZ";
		
	dphi_dz_str[DPHI] = "sdphi";
	dphi_dz_str[DZ] = "sdz";
      }
	
      string meanpointsfile = toad_loader_mean->location(Form("pointsEMC/Points_%s_%s_Mean.txt",Stage_str[istage].c_str(),dPhi_dZ_str[idphi_dz].c_str()) );
      ifstream myfileMeanPoints;
      cout << "Attempting to load points file: " << meanpointsfile.c_str() << endl;
      myfileMeanPoints.open(meanpointsfile.c_str());
	
      string sigmapointsfile = toad_loader_sigma->location( Form("pointsEMC/Points_%s_%s_Sigma.txt",Stage_str[istage].c_str(),dPhi_dZ_str[idphi_dz].c_str()) );
      ifstream myfileSigmaPoints;
      cout << "Attempting to load points file: " << sigmapointsfile.c_str() << endl;
      myfileSigmaPoints.open( sigmapointsfile.c_str()  );

      if (!myfileSigmaPoints)
	{
	  cout << "ERROR!!! points file could not be located!" << endl <<endl;
	}
      else
	{
	  for (int icharge = 0; icharge < NCHARGE; icharge++)
	    {
	      for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++)
		{
		  for (int ized = 0; ized < NZED; ized++)
		    {
		      for (int icent = 0; icent < NCENT; icent++)
			{
									
			  for (int ipt = 0; ipt < NPT; ipt++)
			    {
			      myfileMeanPoints >> mean_points_emc[icharge][iarmsect][ized][icent][ipt][istage][idphi_dz];
			      myfileSigmaPoints >> sigma_points_emc[icharge][iarmsect][ized][icent][ipt][istage][idphi_dz];
			  
			    }
		     
			}//end icent
		    }//end ized
		}//end iarmsect
	    }//end icharge

	  cout << "Successfully loaded points from file" << endl << endl;

	}//end else
	

      myfileMeanPoints.close();
      myfileSigmaPoints.close();

    }//end idphi_dz
  }//end istage  

  delete toad_loader_mean;
  delete toad_loader_sigma;

}//LoadPoints()

void LoadMean()
{

  //initilize the points
 for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++){
   for (int ipt = 0; ipt < NPT; ipt++){
     mean0[iarmsect][ipt] = 0;  
     mean1[iarmsect][ipt] = 0;
     mean2[iarmsect][ipt] = 0;
   }//end ipt
 }//end iarmsect
 
  TOAD *toad_loader = new TOAD("Run13Jet");
  string filename = toad_loader->location("offsetEMC/mean.txt");
  
  ifstream fmean;
  fmean.open(filename.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << filename.c_str() << " could not be located!" << endl <<endl;
    }
  else{
    cout << "Successfully loaded mean from  file " << filename.c_str() << endl << endl;
  }
  for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
	{
          fmean >> mean0[iarmsect][ipt]
		>> mean1[iarmsect][ipt]
		>> mean2[iarmsect][ipt];
	}
    }

  fmean.close();
  delete toad_loader;
  
}//end loadmean

void LoadSigma(){

  //initialize the variables
  for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++){
    for (int ipt = 0; ipt < NPT; ipt++){
      sigma0[iarmsect][ipt] = 0;  
      sigma1[iarmsect][ipt] = 0;
      sigma2[iarmsect][ipt] = 0;
    }//end ipt
  }//end iarmsect
  
  TOAD *toad_loader2 = new TOAD("Run13Jet");
  string filename = toad_loader2->location("offsetEMC/sigma.txt");
  
  cout << "Attempting to load sigma from file " << filename << endl;

  ifstream fsigma;
  fsigma.open(filename.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << filename.c_str() << " could not be located!" << endl <<endl;
    }
  else{
    cout << "Successfully loaded sigma from  file " << filename.c_str() << endl << endl;
  }
  for (int iarmsect = 0; iarmsect < NARMSECTS; iarmsect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          fsigma >> sigma0[iarmsect][ipt]
		 >> sigma1[iarmsect][ipt]
		 >> sigma2[iarmsect][ipt];
        }
    }

  fsigma.close();
  delete toad_loader2;

}//LoadSigma()

float CalcMean(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ){

  if(pt>=2.98) pt = 2.98;
  int ipt = GetPtBin(pt);
  float pt_bin_center = GetPtBinCenter(ipt);
  	
  float mean = mean_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
  /*
  if(ipt>NPT-N_HIGH_PT_TEST_BINS){
    cout << "EmcMatching CalcMean Interp, incorrect ipt, ipt: " << ipt << ", pt: " << pt << endl;
  }
  */	
  if((ipt==0 || pt>=pt_bin_center) && ipt!=NPT-N_HIGH_PT_TEST_BINS-1){
    float y0 = mean_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
    float y1 = mean_points_emc[icharge][iarmsect][ized][icent][ipt+1][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenter(ipt);
    float x1 = GetPtBinCenter(ipt+1);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
    
    mean = (pt - x0)*slope+y0;
  }
  else if(ipt==NPT-N_HIGH_PT_TEST_BINS-1 || pt<pt_bin_center){
    float y0 = mean_points_emc[icharge][iarmsect][ized][icent][ipt-1][STAGE][DPHI_DZ];
    float y1 = mean_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenter(ipt-1);
    float x1 = GetPtBinCenter(ipt);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    mean = (pt - x0)*slope+y0;
  }
  else{
    cout << "EmcMatching CalcMean Interp Broken, ipt: " << ipt << endl;
  }

  return mean;

}//end calcmean()

float CalcSigma(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ){

  if(pt>=2.98) pt = 2.98;
  int ipt = GetPtBin(pt);
  float pt_bin_center = GetPtBinCenter(ipt);
  //  float pt_bin_center_plus1 = GetPtBinCenter(ipt+1);
	
  float sigma = sigma_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
  /*
  if(ipt>NPT-N_HIGH_PT_TEST_BINS){
    cout << "EmcMatching CalcSigma Interp, incorrect ipt, ipt: " << ipt << ", pt: " << pt << endl;
  }
  */	
  if((ipt==0 || pt>=pt_bin_center) && ipt!=NPT-N_HIGH_PT_TEST_BINS-1){
    float y0 = sigma_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
    float y1 = sigma_points_emc[icharge][iarmsect][ized][icent][ipt+1][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenter(ipt);
    float x1 = GetPtBinCenter(ipt+1);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    sigma = (pt - x0)*slope + y0;
  }
  else if(ipt==NPT-N_HIGH_PT_TEST_BINS-1 || pt<pt_bin_center){
    float y0 = sigma_points_emc[icharge][iarmsect][ized][icent][ipt-1][STAGE][DPHI_DZ];
    float y1 = sigma_points_emc[icharge][iarmsect][ized][icent][ipt][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenter(ipt-1);
    float x1 = GetPtBinCenter(ipt);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    sigma = (pt - x0)*slope+y0;
  }
  else{
    cout << "EmcMatching CalcSigma Interp Broken, ipt: " << ipt << endl;
  }
	
  return sigma;

}//end calcsigma()

float CalcEMC(float pt, int icharge, int iarmsect, int ized, int icent, float dphi_dz, int STAGE, int DPHI_DZ)
{
  if (pt > 3.0) pt = 3.0;

  float mean = CalcMean(pt, icharge, iarmsect, ized, icent, STAGE, DPHI_DZ);
  float sigma = CalcSigma(pt, icharge, iarmsect, ized, icent, STAGE, DPHI_DZ);
  
  float sigmalized = (dphi_dz - mean) / sigma;
  
  return sigmalized;

}//end calc()
