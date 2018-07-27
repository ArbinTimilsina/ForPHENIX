#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>

using namespace std;

#define NARMSECT 8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72
#define NECOREBIN 8
#define NSIGMA 2

int IsPbGl(int armsect)
{
  return ((armsect == 4 || armsect == 5) ? 1 : 0);
}


int IsValidYZ(int as, int y, int z)
{
  int ret = 0;
  if (IsPbGl(as))
    {
      if (y >= 0 && y < YPOS_PBGL && z >= 0 && z < ZPOS_PBGL) ret = 1;
    }
  else
    {
      if (y >= 0 && y < YPOS_PBSC && z >= 0 && z < ZPOS_PBSC)   ret = 1;
    }
  return ret;
}

void getStatistics_hot()
{
  int start = 3;
  gStyle->SetOptStat(0);

  const char *armsect[8] = {"Sector: 0(PbSc)", "Sector: 1(PbSc)", "Sector: 2(PbSc)", "Sector: 3(PbSc)",
                            "Sector: 4(PbGl)", "Sector: 5(PbGl)", "Sector: 6(PbSc)", "Sector: 7(PbSc)"
  };

  const char *ebin[8] = {"0.5", "0.7", "1.0", "1.5", "2.0", "3.0", "5.0", "7.0"};

  TFile *file = new TFile("/direct/phenix+hhj/arbint/RootFiles/CuAuEMCal.root", "r");

  //Get the 3D histogram
  TH3D *hEcore[NARMSECT];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      char hname[256];
      sprintf(hname, "hEmc3D_Armsect_%i", ias);
      hEcore[ias] = (TH3D*)file->Get(hname);
    }

  //Read in mean and sigma
  float mean[NARMSECT][NECOREBIN];
  float sigma[NARMSECT][NECOREBIN];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          mean[ias][iec]  = 0.0;
	  sigma[ias][iec] = 0.0;
        }
    }

  cout << "Reading the mean and sigma " << endl;
  ifstream fmean;
  fmean.open("meanSigma.txt");

  if (!fmean)
    {
      cout << "ERROR!!! File could not be located!" << endl;
    }


  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          fmean >> mean[ias][iec]
		>> sigma[ias][iec];
        }
    }

  fmean.close();

  int total[NARMSECT][NECOREBIN][NSIGMA];
  int hot[NARMSECT][NECOREBIN][NSIGMA];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
	  for (int isg = 0; isg < NSIGMA; isg++)
	    {
	      total[ias][iec][isg] = 0;
	      hot[ias][iec][isg] = 0;
	    }
	}
    }

  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          for (int iz = 0; iz < ZPOS_PBGL; iz++ )
            {
              for (int iy = 0; iy < YPOS_PBGL; iy++ )
                {
                  if (IsValidYZ(ias, iy, iz))
                    {
                      float nhitHot = 0.0;
                      for (int iecc = iec; iecc < NECOREBIN; iecc++)
                        {
                          nhitHot += hEcore[ias]->GetBinContent(hEcore[ias]->GetXaxis()->FindBin(iz),
                                                                hEcore[ias]->GetYaxis()->FindBin(iy),
                                                                iecc+1);
                        }

                      for (int isg = 0; isg < NSIGMA; isg++)
                        {
                          float limit_high = mean[ias][iec] + (isg + start) * sigma[ias][iec];
                          total[ias][iec][isg]++;
                          if (nhitHot > limit_high) hot[ias][iec][isg]++;
                        }
                    }
                }
            }
        }
    }

  TH1D *hStatistics[NARMSECT][NSIGMA];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int isg = 0; isg < NSIGMA; isg++)
        {
          char hname[256];
          sprintf(hname, "hStatistics_Armsect%i_sigma%i", ias, isg);

          hStatistics[ias][isg] = new TH1D(hname, hname, NECOREBIN, 0, NECOREBIN);
          hStatistics[ias][isg]->SetTitle(Form("%s", armsect[ias]));
          hStatistics[ias][isg]->SetMarkerStyle(28);
          hStatistics[ias][isg]->SetMarkerColor(isg + 1);
          if (isg == 4) hStatistics[ias][isg]->SetMarkerColor(6);
	  hStatistics[ias][isg]->GetYaxis()->SetRangeUser(0, 10);
          hStatistics[ias][isg]->GetYaxis()->SetTitle("% of hot towers");
          hStatistics[ias][isg]->GetXaxis()->SetLabelSize(0.05);
          hStatistics[ias][isg]->GetYaxis()->SetTitleSize(0.035);

          for (int iec = 0; iec < NECOREBIN; iec++)
            {
              hStatistics[ias][isg]->GetXaxis()->SetBinLabel(iec + 1, ebin[iec]);

	      hStatistics[ias][isg]->SetBinContent(iec + 1,(float)(hot[ias][iec][isg])*100/total[ias][iec][isg]);
            }
        }
    }

  //Plot the Statistics histogram
  TCanvas *cstats = new TCanvas("cstats", "cstats", 1250, 700);
  TLegend *legend[NARMSECT];
  cstats->Divide(4, 2);
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      legend[ias] = new TLegend(0.65, 0.71, 0.85, 0.88);
      legend[ias]->SetBorderSize(0);
      legend[ias]->SetLineColor(1);
      legend[ias]->SetLineStyle(1);
      legend[ias]->SetLineWidth(1);
      legend[ias]->SetFillColor(10);
      legend[ias]->SetTextFont(42);
      legend[ias]->SetTextSize(0.038);

      cstats->cd(ias + 1);
      for (int isg = 0; isg < NSIGMA; isg++)
        {
          if (isg == 0)
            {
              hStatistics[ias][isg]->Draw("p");
            }
          else
            {
              hStatistics[ias][isg]->Draw("pSAME");
            }

          legend[ias]->AddEntry(hStatistics[ias][isg], Form("%i #sigma cut", (isg + start)), "p");
          legend[ias]->Draw();
        }
    }

  cstats->SaveAs("/direct/phenix+hhj/arbint/plots/WarnMaps/statistics_hot.png");
  //  cstats->Close();
}

