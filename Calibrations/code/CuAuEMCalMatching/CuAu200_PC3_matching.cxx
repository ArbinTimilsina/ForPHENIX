#include "CuAu200_PC3_matching.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "TOAD.h"

using namespace std;

void pc3_init_fit_pars_I()
{
  //READ PARAMETERS OUT TO A TEXT FILE
  char inTitle[200];
  TOAD *toad_loader1 = new TOAD("forEMCalMatching");
  for (int iCent = 0; iCent < pc3_cent_bins; iCent++)
    {
      sprintf(inTitle, "%s%i%s", "parametersPC3/CuAu200_pc3_par_C", iCent, "_I.txt");
      string file_location1 = toad_loader1->location(inTitle);
      
      ifstream infile(file_location1.c_str());
      if (!infile)
        {
          cout << "Could not open input PC3 Matching Parameter file " << inTitle << " !!" << endl;
	  delete toad_loader1;
          return ;
        }
      for (int iDet = 0; iDet < pc3_num_det; iDet++)
        {
          for (int iZed = 0; iZed < pc3_zed_bins; iZed++)
            {
              for (int i = 0; i < pc3_func_deg_mean_I; i++)
                {
                  infile >> pc3_zed_mean_I[iCent][iDet][iZed][i];
                }
              for (int i = 0; i < pc3_func_deg_sigma_I; i++)
                {
                  infile >> pc3_zed_sigma_I[iCent][iDet][iZed][i];
                }
            }
        }
      infile.close();
    }
  delete toad_loader1;
  return;
}


void pc3_init_fit_pars_II()
{
  pc3_init_fit_pars_I();
  //READ PARAMETERS OUT TO A TEXT FILE
  char inTitle[200];
  TOAD *toad_loader2 = new TOAD("forEMCalMatching");
  for (int iCent = 0; iCent < pc3_cent_bins; iCent++)
    {
      sprintf(inTitle, "%s%i%s", "parametersPC3/CuAu200_pc3_par_C", iCent, "_II.txt");
      string file_location2 = toad_loader2->location(inTitle);

      ifstream infile(file_location2.c_str());
      if (!infile)
        {
          cout << "Could not open input PC3 Matching Parameter II file " << inTitle << " !!" << endl;
	  delete toad_loader2;
          return ;
        }
      for (int iDet = 0; iDet < pc3_num_det; iDet++)
        {
          for (int iZed = 0; iZed < pc3_zed_bins; iZed++)
            {
              for (int i = 0; i < pc3_func_deg_mean_II; i++)
                {
                  infile >> pc3_zed_mean_II[iCent][iDet][iZed][i];
                }
              for (int i = 0; i < pc3_func_deg_sigma_II; i++)
                {
                  infile >> pc3_zed_sigma_II[iCent][iDet][iZed][i];
                }
            }
        }
      infile.close();
    }
  delete toad_loader2;
  return;
}


float pc3_sdphi_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi)
{
  if (iCent > 4)    iCent = 4.0;
  if (pt  > 3.5)      pt  = 3.5;
  else if (pt < 0.5) pt   = 0.5;

  if     (pc3dphi <= -9999)         return -9999;
  else if (dcarm == 0 && charge > 0)   return (pc3dphi - pc3_eval_mean_I(pt, 0, iCent, iZed)) / pc3_eval_sigma_I(pt, 0, iCent, iZed);
  else if (dcarm == 0 && charge < 0)   return (pc3dphi - pc3_eval_mean_I(pt, 4, iCent, iZed)) / pc3_eval_sigma_I(pt, 4, iCent, iZed);
  else if (dcarm == 1 && charge > 0)   return (pc3dphi - pc3_eval_mean_I(pt, 1, iCent, iZed)) / pc3_eval_sigma_I(pt, 1, iCent, iZed);
  else if (dcarm == 1 && charge < 0)   return (pc3dphi - pc3_eval_mean_I(pt, 5, iCent, iZed)) / pc3_eval_sigma_I(pt, 5, iCent, iZed);
  else                            return -9999;
}


float pc3_sdz_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz)
{
  if (iCent > 4)    iCent = 4.0;
  if (pt > 3.5)      pt   = 3.5;
  else if (pt < 0.5) pt   = 0.5;
  if (charge < 0  &&  dcarm == 0  &&  iCent >= 4)  iCent = 3;   //special condition for poor statistics

  if     (pc3dz <= -9999)          return -9999;
  else if (dcarm == 0 && charge > 0)  return (pc3dz - pc3_eval_mean_I(pt, 2, iCent, iZed)) / pc3_eval_sigma_I(pt, 2, iCent, iZed);
  else if (dcarm == 0 && charge < 0)  return (pc3dz - pc3_eval_mean_I(pt, 6, iCent, iZed)) / pc3_eval_sigma_I(pt, 6, iCent, iZed);
  else if (dcarm == 1 && charge > 0)  return (pc3dz - pc3_eval_mean_I(pt, 3, iCent, iZed)) / pc3_eval_sigma_I(pt, 3, iCent, iZed);
  else if (dcarm == 1 && charge < 0)  return (pc3dz - pc3_eval_mean_I(pt, 7, iCent, iZed)) / pc3_eval_sigma_I(pt, 7, iCent, iZed);
  else                           return -9999;
}


float pc3_sdphi_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi)
{
  float pc3sdphi = pc3_sdphi_func_I(charge, dcarm, iCent, iZed, pt, pc3dphi);

  if (iCent > 4)    iCent = 4;
  if (pt > 2.0)      pt = 2.0;
  //  else if(pt < 0.5) pt = 0.5;

  if     (dcarm == 0 && charge > 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 0, iCent, iZed)) / pc3_eval_sigma_II(pt, 0, iCent, iZed);
  else if (dcarm == 0 && charge < 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 4, iCent, iZed)) / pc3_eval_sigma_II(pt, 4, iCent, iZed);
  else if (dcarm == 1 && charge > 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 1, iCent, iZed)) / pc3_eval_sigma_II(pt, 1, iCent, iZed);
  else if (dcarm == 1 && charge < 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 5, iCent, iZed)) / pc3_eval_sigma_II(pt, 5, iCent, iZed);
  else return -9999;
}

float pc3_sdz_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz)
{
  float pc3sdz = pc3_sdz_func_I(charge, dcarm, iCent, iZed, pt, pc3dz);

  if (iCent > 4)    iCent = 4;
  if (pt > 2.0)      pt = 2.0;
  //  else if(pt < 0.5) pt = 0.5;

  if     (dcarm == 0 && charge > 0)  return (pc3sdz - pc3_eval_mean_II(pt, 2, iCent, iZed)) / pc3_eval_sigma_II(pt, 2, iCent, iZed);
  else if (dcarm == 0 && charge < 0)  return (pc3sdz - pc3_eval_mean_II(pt, 6, iCent, iZed)) / pc3_eval_sigma_II(pt, 6, iCent, iZed);
  else if (dcarm == 1 && charge > 0)  return (pc3sdz - pc3_eval_mean_II(pt, 3, iCent, iZed)) / pc3_eval_sigma_II(pt, 3, iCent, iZed);
  else if (dcarm == 1 && charge < 0)  return (pc3sdz - pc3_eval_mean_II(pt, 7, iCent, iZed)) / pc3_eval_sigma_II(pt, 7, iCent, iZed);
  else return -9999;
}




float pc3_eval_mean_I(float pt, int iDet, int iCent, int iZed)
{
  float mean = pc3_zed_mean_I[iCent][iDet][iZed][0]            + pc3_zed_mean_I[iCent][iDet][iZed][1] * pt
    + pc3_zed_mean_I[iCent][iDet][iZed][2] * pow(pt, 2)  + pc3_zed_mean_I[iCent][iDet][iZed][3] * pow(pt, 3)
    + pc3_zed_mean_I[iCent][iDet][iZed][4] * pow(pt, 4)  + pc3_zed_mean_I[iCent][iDet][iZed][5] * pow(pt, 5)
    + pc3_zed_mean_I[iCent][iDet][iZed][6] / sqrt(pt);
  return mean;
}

float pc3_eval_sigma_I(float pt, int iDet, int iCent, int iZed)
{
  float sigma = pc3_zed_sigma_I[iCent][iDet][iZed][0]           + pc3_zed_sigma_I[iCent][iDet][iZed][1] * pt
    + pc3_zed_sigma_I[iCent][iDet][iZed][2] / pt        + pc3_zed_sigma_I[iCent][iDet][iZed][3] * pow(pt, 2)
    + pc3_zed_sigma_I[iCent][iDet][iZed][4] * pow(pt, 3) + pc3_zed_sigma_I[iCent][iDet][iZed][5] * pow(pt, 4)
    + pc3_zed_sigma_I[iCent][iDet][iZed][6] * pow(pt, 5) + pc3_zed_sigma_I[iCent][iDet][iZed][7] / sqrt(pt);
  return sigma;
}

float pc3_eval_mean_II(float pt, int iDet, int iCent, int iZed)
{
  float mean = pc3_zed_mean_II[iCent][iDet][iZed][0]          + pc3_zed_mean_II[iCent][iDet][iZed][1] * pt
               + pc3_zed_mean_II[iCent][iDet][iZed][2] / pt       + pc3_zed_mean_II[iCent][iDet][iZed][3] * pow(pt, 2)
               + pc3_zed_mean_II[iCent][iDet][iZed][4] / sqrt(pt);

  return mean;
}

float pc3_eval_sigma_II(float pt, int iDet, int iCent, int iZed)
{
  float sigma = pc3_zed_sigma_II[iCent][iDet][iZed][0]           + pc3_zed_sigma_II[iCent][iDet][iZed][1] * pt
                + pc3_zed_sigma_II[iCent][iDet][iZed][2] / pt        + pc3_zed_sigma_II[iCent][iDet][iZed][3] * pow(pt, 2)
                + pc3_zed_sigma_II[iCent][iDet][iZed][4] * pow(pt, 3) + pc3_zed_sigma_II[iCent][iDet][iZed][5] * pow(pt, 4)
                + pc3_zed_sigma_II[iCent][iDet][iZed][6] / sqrt(pt);
  return sigma;
}
