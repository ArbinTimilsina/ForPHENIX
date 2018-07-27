//#include <CommonVar.h>
#include <PC3Matching.h>
#include <iostream>
#include <string>
#include <fstream>

#include <TOAD.h>
#include <TString.h>

using namespace std;

enum type1_enum{DPHI, DZ};
enum type2_enum{INITIAL, INTERMEDIATE, FINAL};

const double  cuthigh[NPT_PC3] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 10.0, 25.0};
const double  cutlow[NPT_PC3] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 10.0};

int GetPtBinPC3(float pT){
  int kbin = -1;
  for (int m1 = 0; m1 < NPT_PC3; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

float GetPtBinCenterPC3(int pt_bin)
{
  float pt_bin_center = (cutlow[pt_bin]+cuthigh[pt_bin])/2.;
  return pt_bin_center;
}

void LoadAllPC3Files(){

  //for sdphi and sdz
  LoadPointsPC3();

}//end LoadAllFiles()

void LoadPointsPC3(){

  //initialize the points
  for (int icharge = 0; icharge < NCHARGE_PC3; icharge++){
    for (int iarm = 0; iarm < NARMS_PC3; iarm++){
      for (int ized = 0; ized < NZED_PC3; ized++){
	for (int icent = 0; icent < NCENT_PC3; icent++){
	  for (int ipt = 0; ipt < NPT_PC3; ipt++){
	    for(int istage = 0; istage < 3; istage++){
	      for(int idphi_dz = 0; idphi_dz < 2; idphi_dz++){
		
		mean_points_pc3[icharge][iarm][ized][icent][ipt][istage][idphi_dz] = -9999.0;
		sigma_points_pc3[icharge][iarm][ized][icent][ipt][istage][idphi_dz] = -9999.0;
		
	      }//end idphi_dz
	    }//end istage    
	  }//end ipt
	}
      }
    }
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
	
      string meanpointsfile = toad_loader_mean->location(Form("pointsPC3/Points_%s_%s_Mean.txt",Stage_str[istage].c_str(),dPhi_dZ_str[idphi_dz].c_str()) );
      ifstream myfileMeanPoints;
      cout << "Attempting to load points file: " << meanpointsfile.c_str() << endl;
      myfileMeanPoints.open(meanpointsfile.c_str());
	
      string sigmapointsfile = toad_loader_sigma->location( Form("pointsPC3/Points_%s_%s_Sigma.txt",Stage_str[istage].c_str(),dPhi_dZ_str[idphi_dz].c_str()) );
      ifstream myfileSigmaPoints;
      cout << "Attempting to load points file: " << sigmapointsfile.c_str() << endl;
      myfileSigmaPoints.open( sigmapointsfile.c_str()  );

      if (!myfileSigmaPoints)
	{
	  cout << "ERROR!!! points file could not be located!" << endl <<endl;
	}
      else
	{
	  for (int icharge = 0; icharge < NCHARGE_PC3; icharge++)
	    {
	      for (int iarm = 0; iarm < NARMS_PC3; iarm++)
		{
		  for (int ized = 0; ized < NZED_PC3; ized++)
		    {
		      for (int icent = 0; icent < NCENT_PC3; icent++)
			{
									
			  for (int ipt = 0; ipt < NPT_PC3; ipt++)
			    {
			      myfileMeanPoints >> mean_points_pc3[icharge][iarm][ized][icent][ipt][istage][idphi_dz];
			      myfileSigmaPoints >> sigma_points_pc3[icharge][iarm][ized][icent][ipt][istage][idphi_dz];
			  
			    }
		     
			}//end icent
		    }//end ized
		}//end iarm
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

float CalcMeanPC3(float pt, int icharge, int iarm, int ized, int icent, int STAGE, int DPHI_DZ){

  if(pt>=2.98) pt = 2.98;
  int ipt = GetPtBinPC3(pt);
  float pt_bin_center = GetPtBinCenterPC3(ipt);
  	
  float mean = mean_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
  /*
  if(ipt>=NPT_PC3-N_HIGH_PT_TEST_BINS_PC3){
    cout << "PC3Matching CalcMean Interp, incorrect ipt, ipt: " << ipt << ", pt: " << pt << endl;
  }
  */	
  if((ipt==0 || pt>=pt_bin_center) && ipt!=NPT_PC3-N_HIGH_PT_TEST_BINS_PC3-1){
    float y0 = mean_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
    float y1 = mean_points_pc3[icharge][iarm][ized][icent][ipt+1][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenterPC3(ipt);
    float x1 = GetPtBinCenterPC3(ipt+1);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
    
    mean = (pt - x0)*slope+y0;
  }
  else if(ipt==NPT_PC3-N_HIGH_PT_TEST_BINS_PC3-1 || pt<pt_bin_center){
    float y0 = mean_points_pc3[icharge][iarm][ized][icent][ipt-1][STAGE][DPHI_DZ];
    float y1 = mean_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenterPC3(ipt-1);
    float x1 = GetPtBinCenterPC3(ipt);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    mean = (pt - x0)*slope+y0;
  }
  else{
    cout << "PC3Matching CalcMean Interp Broken, ipt: " << ipt << endl;
  }

  return mean;

}//end calcmean()

float CalcSigmaPC3(float pt, int icharge, int iarm, int ized, int icent, int STAGE, int DPHI_DZ){

  if(pt>=2.98) pt = 2.98;
  int ipt = GetPtBinPC3(pt);
  float pt_bin_center = GetPtBinCenterPC3(ipt);
  //  float pt_bin_center_plus1 = GetPtBinCenterPC3(ipt+1);
	
  float sigma = sigma_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
  /*
  if(ipt>=NPT_PC3-N_HIGH_PT_TEST_BINS_PC3){
    cout << "PC3Matching CalcSigma Interp, incorrect ipt, ipt: " << ipt << ", pt: " << pt << endl;
  }
  */	
  if((ipt==0 || pt>=pt_bin_center) && ipt!=NPT_PC3-N_HIGH_PT_TEST_BINS_PC3-1){
    float y0 = sigma_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
    float y1 = sigma_points_pc3[icharge][iarm][ized][icent][ipt+1][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenterPC3(ipt);
    float x1 = GetPtBinCenterPC3(ipt+1);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    sigma = (pt - x0)*slope + y0;
  }
  else if(ipt==NPT_PC3-N_HIGH_PT_TEST_BINS_PC3-1 || pt<pt_bin_center){
    float y0 = sigma_points_pc3[icharge][iarm][ized][icent][ipt-1][STAGE][DPHI_DZ];
    float y1 = sigma_points_pc3[icharge][iarm][ized][icent][ipt][STAGE][DPHI_DZ];
		
    float x0 = GetPtBinCenterPC3(ipt-1);
    float x1 = GetPtBinCenterPC3(ipt);
		
    float slope_numer = y1 - y0;
    float slope_denom = x1 - x0;
    float slope = slope_numer/slope_denom;
		
    sigma = (pt - x0)*slope+y0;
  }
  else{
    cout << "PC3Matching CalcSigma Interp Broken, ipt: " << ipt << endl;
  }
	
  return sigma;

}//end calcsigma()

float CalcPC3(float pt, int icharge, int iarm, int ized, int icent, float dphi_dz, int STAGE, int DPHI_DZ)
{
  if (pt > 3.0) pt = 3.0;

  float mean = CalcMeanPC3(pt, icharge, iarm, ized, icent, STAGE, DPHI_DZ);
  float sigma = CalcSigmaPC3(pt, icharge, iarm, ized, icent, STAGE, DPHI_DZ);
  
  float sigmalized = (dphi_dz - mean) / sigma;
  
  return sigmalized;

}//end calc()
