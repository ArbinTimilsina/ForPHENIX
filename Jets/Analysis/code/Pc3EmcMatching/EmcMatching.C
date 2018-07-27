#include "EmcMatching.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "TOAD.h"

using namespace std;

double cutlow[NPT] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0};
double cuthigh[NPT] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};

int GetPtBin(float pT)
{
  int kbin = -1;
  for (int m1 = 0; m1 < NPT; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

void LoadEmcMatchingParameters()
{
  LoadInitialdPhiMeanParameters();
  LoadInitialdPhiSigmaParameters();

  LoadFinaldPhiMeanParameters();
  LoadFinaldPhiSigmaParameters();

  LoadMean();
  LoadSigma();

  LoadInitialdZMeanParameters();
  LoadInitialdZSigmaParameters();

  LoadFinaldZMeanParameters();
  LoadFinaldZSigmaParameters();
}


void LoadInitialdPhiMeanParameters()
{
  TOAD *toad_loader1 = new TOAD("Pc3EmcMatching");
  string initialdPhiMean = toad_loader1->location("parametersEMC/Parameters_Initial_dPhi_Mean.txt");
  cout << "Loading initial dphi mean parameters from file " << initialdPhiMean << endl;

  ifstream fmean;
  fmean.open(initialdPhiMean.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << initialdPhiMean << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fmean >> initialdPhiMeanParam_0[icharge][isect][ized][icent]
			>> initialdPhiMeanParam_1[icharge][isect][ized][icent]
			>> initialdPhiMeanParam_2[icharge][isect][ized][icent];
                }
            }
        }
    }
  fmean.close();
  delete toad_loader1;
  cout << "Successfully loaded initial dphi mean parameters from  file " << initialdPhiMean << endl << endl;
}


void LoadInitialdPhiSigmaParameters()
{
  TOAD *toad_loader2 = new TOAD("Pc3EmcMatching");
  string initialdPhiSigma = toad_loader2->location("parametersEMC/Parameters_Initial_dPhi_Sigma.txt");
  cout << "Loading initial dphi sigma parameters from file " << initialdPhiSigma << endl;

  ifstream fsigma;
  fsigma.open(initialdPhiSigma.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << initialdPhiSigma << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fsigma >> initialdPhiSigmaParam_0[icharge][isect][ized][icent]
			 >> initialdPhiSigmaParam_1[icharge][isect][ized][icent]
			 >> initialdPhiSigmaParam_2[icharge][isect][ized][icent];
                }
            }
        }
    }
  fsigma.close();
  delete toad_loader2;
  cout << "Successfully loaded initial dphi sigma parameters from  file " << initialdPhiSigma << endl << endl;
}

void LoadFinaldPhiMeanParameters()
{
  TOAD *toad_loader3 = new TOAD("Pc3EmcMatching");
  string finaldPhiMean = toad_loader3->location("parametersEMC/Parameters_Final_dPhi_Mean.txt");
  cout << "Loading final dphi mean parameters from file " << finaldPhiMean << endl;

  ifstream fmean;
  fmean.open(finaldPhiMean.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << finaldPhiMean << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fmean >> finaldPhiMeanParam_0[icharge][isect][ized][icent]
			>> finaldPhiMeanParam_1[icharge][isect][ized][icent]
			>> finaldPhiMeanParam_2[icharge][isect][ized][icent];
                }
            }
        }
    }
  fmean.close();
  delete toad_loader3;
  cout << "Successfully loaded final dphi mean parameters from  file " << finaldPhiMean << endl << endl;
}

void LoadFinaldPhiSigmaParameters()
{
  TOAD *toad_loader4 = new TOAD("Pc3EmcMatching");
  string finaldPhiSigma = toad_loader4->location("parametersEMC/Parameters_Final_dPhi_Sigma.txt");
  cout << "Loading final dphi sigma parameters from file " << finaldPhiSigma << endl;

  ifstream fsigma;
  fsigma.open(finaldPhiSigma.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << finaldPhiSigma << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fsigma >> finaldPhiSigmaParam_0[icharge][isect][ized][icent];
                }
            }
        }
    }
  fsigma.close();
  delete toad_loader4;
  cout << "Successfully loaded final dphi sigma parameters from  file " << finaldPhiSigma << endl << endl;
}

float CalculateInitialEmcsdPhi(int charge, int sector, int zed, int cent, float emcdphi, float pt)
{
  float mean_p0 = initialdPhiMeanParam_0[charge][sector][zed][cent];
  float mean_p1 = initialdPhiMeanParam_1[charge][sector][zed][cent];
  float mean_p2 = initialdPhiMeanParam_2[charge][sector][zed][cent];

  float sigma_p0 = initialdPhiSigmaParam_0[charge][sector][zed][cent];
  float sigma_p1 = initialdPhiSigmaParam_1[charge][sector][zed][cent];
  float sigma_p2 = initialdPhiSigmaParam_2[charge][sector][zed][cent];

  if (pt > 3.0) pt = 3.0;

  float mean = mean_p0 + mean_p1 / pt + mean_p2 / pt / pt;
  float sigma = sigma_p0 + sigma_p1 / pt + sigma_p2 / pt / pt;

  float emcsdphi = (emcdphi - mean) / sigma;

  return emcsdphi;
}


float CalculateFinalEmcsdPhi(int charge, int sector, int zed, int cent, float emcsdphi, float pt)
{
  float mean_p0 = finaldPhiMeanParam_0[charge][sector][zed][cent];
  float mean_p1 = finaldPhiMeanParam_1[charge][sector][zed][cent];
  float mean_p2 = finaldPhiMeanParam_2[charge][sector][zed][cent];

  float sigma_p0 = finaldPhiSigmaParam_0[charge][sector][zed][cent];

  if (pt > 3.0) pt = 3.0;

  float mean = mean_p0 + mean_p1 / pt + mean_p2 / pt / pt;
  float sigma = sigma_p0;

  float emcsdphi_final = (emcsdphi - mean) / sigma;

  return emcsdphi_final;
}



void LoadMean()
{
  TOAD *toad_loader5 = new TOAD("Pc3EmcMatching");
  string mean = toad_loader5->location("offset/mean.txt");
  cout << "Loading mean parameters from file " << mean << endl;

  ifstream fmean;
  fmean.open(mean.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << mean << " could not be located!" << endl;
    }

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          fmean >> mean0[isect][ipt]
                >> mean1[isect][ipt]
                >> mean2[isect][ipt];
        }
    }

  fmean.close();
  delete toad_loader5;
  cout << "Successfully loaded mean parameters from  file " << mean << endl << endl;

}

void LoadSigma()
{
  TOAD *toad_loader6 = new TOAD("Pc3EmcMatching");
  string sigma = toad_loader6->location("offset/sigma.txt");
  cout << "Loading sigma parameters from file " << sigma << endl;

  ifstream fsigma;
  fsigma.open(sigma.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << sigma << " could not be located!" << endl;
    }

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
	{
          fsigma >> sigma0[isect][ipt]
                 >> sigma1[isect][ipt]
                 >> sigma2[isect][ipt];
        }
    }

  fsigma.close();
  delete toad_loader6;
  cout << "Successfully loaded sigma parameters from  file " << sigma << endl << endl;


}

void LoadInitialdZMeanParameters()
{
  TOAD *toad_loader5 = new TOAD("Pc3EmcMatching");
  string initialdZMean = toad_loader5->location("parametersEMC/Parameters_Initial_dZ_Mean.txt");
  cout << "Loading initial dz mean parameters from file " << initialdZMean << endl;

  ifstream fmean;
  fmean.open(initialdZMean.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << initialdZMean << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fmean >> initialdZMeanParam_0[icharge][isect][ized][icent];
                }
            }
        }
    }
  fmean.close();
  delete toad_loader5;
  cout << "Successfully loaded initial dz mean parameters from  file " << initialdZMean << endl << endl;
}

void LoadInitialdZSigmaParameters()
{
  TOAD *toad_loader6 = new TOAD("Pc3EmcMatching");
  string initialdZSigma = toad_loader6->location("parametersEMC/Parameters_Initial_dZ_Sigma.txt");
  cout << "Loading initial dz sigma parameters from file " << initialdZSigma << endl;

  ifstream fsigma;
  fsigma.open(initialdZSigma.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << initialdZSigma << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fsigma >> initialdZSigmaParam_0[icharge][isect][ized][icent]
			 >> initialdZSigmaParam_1[icharge][isect][ized][icent]
			 >> initialdZSigmaParam_2[icharge][isect][ized][icent];
                }
            }
        }
    }
  fsigma.close();
  delete toad_loader6;
  cout << "Successfully loaded initial dz sigma parameters from  file " << initialdZSigma << endl << endl;
}

void LoadFinaldZMeanParameters()
{
  TOAD *toad_loader7 = new TOAD("Pc3EmcMatching");
  string finaldZMean = toad_loader7->location("parametersEMC/Parameters_Final_dZ_Mean.txt");
  cout << "Loading final dz mean parameters from file " << finaldZMean << endl;

  ifstream fmean;
  fmean.open(finaldZMean.c_str());

  if (!fmean)
    {
      cout << "ERROR!!! File " << finaldZMean << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fmean >> finaldZMeanParam_0[icharge][isect][ized][icent];
                }
            }
        }
    }
  fmean.close();
  delete toad_loader7;
  cout << "Successfully loaded final dz mean parameters from  file " << finaldZMean << endl << endl;
}

void LoadFinaldZSigmaParameters()
{
  TOAD *toad_loader8 = new TOAD("Pc3EmcMatching");
  string finaldZSigma = toad_loader8->location("parametersEMC/Parameters_Final_dZ_Sigma.txt");
  cout << "Loading final dz sigma parameters from file " << finaldZSigma << endl;

  ifstream fsigma;
  fsigma.open(finaldZSigma.c_str());

  if (!fsigma)
    {
      cout << "ERROR!!! File " << finaldZSigma << " could not be located!" << endl;
    }

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  fsigma >> finaldZSigmaParam_0[icharge][isect][ized][icent]
			 >> finaldZSigmaParam_1[icharge][isect][ized][icent]
			 >> finaldZSigmaParam_2[icharge][isect][ized][icent];
                }
            }
        }
    }
  fsigma.close();
  delete toad_loader8;
  cout << "Successfully loaded final dz sigma parameters from  file " << finaldZSigma << endl << endl;
}

float CalculateCorrectedEmcdZ(float beta, float pt, float emcdz, int sector)
{
  float tanTheta = tan((3.1416 / 2) - beta);
  int ptBin = GetPtBin(pt);

  float meann = mean0[sector][ptBin] + mean1[sector][ptBin]*tanTheta +  mean2[sector][ptBin]*tanTheta*tanTheta*tanTheta;
  float sigmaa = sigma0[sector][ptBin] + sigma0[sector][ptBin]*tanTheta*tanTheta + sigma1[sector][ptBin]*tanTheta*tanTheta*tanTheta*tanTheta;

  float emcdz_corrected = (emcdz - meann)/(sigmaa);

  return emcdz_corrected;
}

float CalculateInitialEmcsdZ(int charge, int sector, int zed, int cent, float emcdz, float pt)
{
  float mean_p0 = initialdZMeanParam_0[charge][sector][zed][cent];

  float sigma_p0 = initialdZSigmaParam_0[charge][sector][zed][cent];
  float sigma_p1 = initialdZSigmaParam_1[charge][sector][zed][cent];
  float sigma_p2 = initialdZSigmaParam_2[charge][sector][zed][cent];

  if (pt > 3.0) pt = 3.0;

  float mean = mean_p0;
  float sigma = sigma_p0 + sigma_p1 / pt + sigma_p2 /pt /pt;

  float emcsdz = (emcdz - mean) / sigma;

  return emcsdz;
}


float CalculateFinalEmcsdZ(int charge, int sector, int zed, int cent, float emcsdz, float pt)
{
  float mean_p0 = finaldZMeanParam_0[charge][sector][zed][cent];

  float sigma_p0 = finaldZSigmaParam_0[charge][sector][zed][cent];
  float sigma_p1 = finaldZSigmaParam_1[charge][sector][zed][cent];
  float sigma_p2 = finaldZSigmaParam_2[charge][sector][zed][cent];

  if (pt > 3.0) pt = 3.0;

  float mean = mean_p0;
  float sigma = sigma_p0 + sigma_p1 / pt + sigma_p2 / pt / pt;

  float emcsdz_final = (emcsdz - mean) / sigma;

  return emcsdz_final;
}
