#include "Pc3Matching.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "TOAD.h"

using namespace std;

void pc3_init_fit_pars_I()
{
  char inTitle1[200];
  TOAD *toad_loader1 = new TOAD("JetAnalyzer");
  for (int iCent = 0; iCent < pc3_cent_bins; iCent++)
    {
      sprintf(inTitle1, "%s%i%s", "parametersPC3/CuAu200_pc3_par_C", iCent, "_I.txt");
      string file_location1 = toad_loader1->location(inTitle1);
      ifstream infile(file_location1.c_str());

      if (!infile)
        {
          cout << "Could not open input PC3 Matching Parameter file " << inTitle1 << " !!" << endl;
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

  char inTitle2[200];
  TOAD *toad_loader2 = new TOAD("JetAnalyzer");
  for (int iCent = 0; iCent < pc3_cent_bins - 1; iCent++)
    {
      sprintf(inTitle2, "%s%i%s", "parametersPC3/CuAu200_pc3_slopes_C", iCent, "_I.txt");
      string file_location2 = toad_loader2->location(inTitle2);
      ifstream outfile2(file_location2.c_str());
      outfile2.precision(9);

      for (int iDet = 0; iDet < pc3_num_det; iDet++)    //8x
        {
          for (int iZed = 0; iZed < pc3_zed_bins; iZed++)  //10x
            {
              for (int iPt = 0; iPt < pc3_pt_bins - 1; iPt++)   //7x
                {
                  outfile2 >> point_lines_I[iCent][iDet][iZed][iPt][0];
                  outfile2 >> point_lines_I[iCent][iDet][iZed][iPt][1];
                  outfile2 >> point_lines_I[iCent][iDet][iZed][iPt][2];
                  outfile2 >> point_lines_I[iCent][iDet][iZed][iPt][3];
                }
            }
        }
      outfile2.close();
    }
  delete toad_loader2;
  return;
}


void pc3_init_fit_pars_II()
{
  pc3_init_fit_pars_I();

  char inTitle3[200];
  TOAD *toad_loader3 = new TOAD("JetAnalyzer");
  for (int iCent = 0; iCent < pc3_cent_bins; iCent++)
    {
      sprintf(inTitle3, "%s%i%s", "parametersPC3/CuAu200_pc3_par_C", iCent, "_II.txt");
      string file_location3 = toad_loader3->location(inTitle3);
      ifstream infile(file_location3.c_str());

      if (!infile)
        {
          cout << "Could not open input PC3 Matching Parameter II file " << inTitle3 << " !!" << endl;
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
  delete toad_loader3;

  return;
}


float pc3_sdphi_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi)
{

  if (iCent > 6)    iCent = 6;
  if (pt > 3.5)     pt = 3.5;

  if (pc3dphi > -9999)
    {

      if (dcarm == 0 && charge > 0)   return (pc3dphi - pc3_eval_mean_I(pt, 0, iCent, iZed)) / pc3_eval_sigma_I(pt, 0, iCent, iZed);
      else if (dcarm == 0 && charge < 0)   return (pc3dphi - pc3_eval_mean_I(pt, 4, iCent, iZed)) / pc3_eval_sigma_I(pt, 4, iCent, iZed);
      else if (dcarm == 1 && charge > 0)   return (pc3dphi - pc3_eval_mean_I(pt, 1, iCent, iZed)) / pc3_eval_sigma_I(pt, 1, iCent, iZed);
      else if (dcarm == 1 && charge < 0)   return (pc3dphi - pc3_eval_mean_I(pt, 5, iCent, iZed)) / pc3_eval_sigma_I(pt, 5, iCent, iZed);
      else                            return -9999;

    }
  else
    {
      return -9999;
    }
}


float pc3_sdz_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz)
{
  if (iCent > 6)    iCent = 6;
  if (pt > 3.5)      pt = 3.5;

  int iPt = -9999;
  if (pt >= 0.4 && pt < 0.7) iPt = 0;
  else if (pt >= 0.7 && pt < 0.9) iPt = 1;
  else if (pt >= 0.9 && pt < 1.1) iPt = 2;
  else if (pt >= 1.1 && pt < 1.3) iPt = 3;
  else if (pt >= 1.3 && pt < 1.5) iPt = 4;
  else if (pt >= 1.5 && pt < 1.7) iPt = 5;
  else if (pt >= 1.7 && pt < 1.9) iPt = 6;
  else if (pt >= 1.9 && pt < 2.1) iPt = 7;
  else if (pt >= 2.1 && pt < 2.3) iPt = 8;
  else if (pt >= 2.3 && pt < 2.5) iPt = 9;
  else if (pt >= 2.5 && pt < 2.8) iPt = 10;
  else if (pt >= 2.8 && pt < 3.2) iPt = 11;
  else if (pt >= 3.2 && pt < 3.8) iPt = 12;
  else if (pt >= 3.8 && pt < 5.0) iPt = 13;


  if (pc3dz > -9999)
    {

      int iDet = -9999;

      if     (dcarm == 0 && charge > 0)
        {
          iDet = 2;
        }
      else if (dcarm == 0 && charge < 0)
        {
          iDet = 6;
        }
      else if (dcarm == 1 && charge > 0)
        {
          iDet = 3;
        }
      else if (dcarm == 1 && charge < 0)
        {
          iDet = 7;
        }
      else return -9999;

      float mean  = point_lines_I[iCent][iDet][iZed][iPt][0] * pt + point_lines_I[iCent][iDet][iZed][iPt][1];
      float sigma = point_lines_I[iCent][iDet][iZed][iPt][2] * pt + point_lines_I[iCent][iDet][iZed][iPt][3];

      return (pc3dz - mean) / sigma;
    }
  else
    {
      return -9999;
    }

}


float pc3_sdphi_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi)
{
  float pc3sdphi = pc3_sdphi_func_I(charge, dcarm, iCent, iZed, pt, pc3dphi);

  if (iCent > 6)    iCent = 6;
  if (pt > 3.0)      pt = 3.0;

  if     (dcarm == 0 && charge > 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 0, iCent, iZed)) / pc3_eval_sigma_II(pt, 0, iCent, iZed);
  else if (dcarm == 0 && charge < 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 4, iCent, iZed)) / pc3_eval_sigma_II(pt, 4, iCent, iZed);
  else if (dcarm == 1 && charge > 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 1, iCent, iZed)) / pc3_eval_sigma_II(pt, 1, iCent, iZed);
  else if (dcarm == 1 && charge < 0)  return (pc3sdphi - pc3_eval_mean_II(pt, 5, iCent, iZed)) / pc3_eval_sigma_II(pt, 5, iCent, iZed);
  else return -9999;
}


float pc3_sdz_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz)
{
  float pc3sdz = pc3_sdz_func_I(charge, dcarm, iCent, iZed, pt, pc3dz);

  if (iCent > 6)    iCent = 6;
  if (pt > 3.0)      pt = 3.0;

  if     (dcarm == 0 && charge > 0)  return (pc3sdz - pc3_eval_mean_II(pt, 2, iCent, iZed)) / pc3_eval_sigma_II(pt, 2, iCent, iZed);
  else if (dcarm == 0 && charge < 0)  return (pc3sdz - pc3_eval_mean_II(pt, 6, iCent, iZed)) / pc3_eval_sigma_II(pt, 6, iCent, iZed);
  else if (dcarm == 1 && charge > 0)  return (pc3sdz - pc3_eval_mean_II(pt, 3, iCent, iZed)) / pc3_eval_sigma_II(pt, 3, iCent, iZed);
  else if (dcarm == 1 && charge < 0)  return (pc3sdz - pc3_eval_mean_II(pt, 7, iCent, iZed)) / pc3_eval_sigma_II(pt, 7, iCent, iZed);
  else return -9999;
}




float pc3_eval_mean_I(float pt, int iDet, int iCent, int iZed)
{
  float mean = pc3_zed_mean_I[iCent][iDet][iZed][0]              + pc3_zed_mean_I[iCent][iDet][iZed][1] * pt
    + pc3_zed_mean_I[iCent][iDet][iZed][2] * pow(pt, 2)    + pc3_zed_mean_I[iCent][iDet][iZed][3] * pow(pt, 3)
    + pc3_zed_mean_I[iCent][iDet][iZed][4] * pow(pt, 4)    + pc3_zed_mean_I[iCent][iDet][iZed][5] * pow(pt, 5)
    + pc3_zed_mean_I[iCent][iDet][iZed][6] / sqrt(pt);
  return mean;
}

float pc3_eval_sigma_I(float pt, int iDet, int iCent, int iZed)
{
  float sigma = pc3_zed_sigma_I[iCent][iDet][iZed][0]            + pc3_zed_sigma_I[iCent][iDet][iZed][1] * pt
    + pc3_zed_sigma_I[iCent][iDet][iZed][2] / pt         + pc3_zed_sigma_I[iCent][iDet][iZed][3] * pow(pt, 2)
    + pc3_zed_sigma_I[iCent][iDet][iZed][4] * pow(pt, 3)  + pc3_zed_sigma_I[iCent][iDet][iZed][5] * pow(pt, 4)
    + pc3_zed_sigma_I[iCent][iDet][iZed][6] * pow(pt, 5)  + pc3_zed_sigma_I[iCent][iDet][iZed][7] / sqrt(pt);
  return sigma;
}

float pc3_eval_mean_II(float pt, int iDet, int iCent, int iZed)
{
  float mean = pc3_zed_mean_II[iCent][iDet][iZed][0]             + pc3_zed_mean_II[iCent][iDet][iZed][1] * pt
    + pc3_zed_mean_II[iCent][iDet][iZed][2] / pt          + pc3_zed_mean_II[iCent][iDet][iZed][3] * pow(pt, 2)
    + pc3_zed_mean_II[iCent][iDet][iZed][4] / sqrt(pt)    + pc3_zed_mean_II[iCent][iDet][iZed][5] * pow(pt, 3)
    + pc3_zed_mean_II[iCent][iDet][iZed][6] * pow(pt, 4);
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
