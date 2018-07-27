#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH3D.h>
#include <TH1F.h>
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

void getSigma()
{
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);

  const char *armsect[8] = {"Sector: 0(PbSc)", "Sector: 1(PbSc)", "Sector: 2(PbSc)", "Sector: 3(PbSc)", 
			    "Sector: 4(PbGl)", "Sector: 5(PbGl)", "Sector: 6(PbSc)", "Sector: 7(PbSc)"
  };

  const char *ecore[8] = {"ecore >0.5 (GeV/c)", "ecore >0.7 (GeV/c)", "ecore >1.0 (GeV/c)", "ecore >1.5 (GeV/c)", 
			  "ecore >2.0 (GeV/c)", "ecore >3.0 (GeV/c)", "ecore >5.0 (GeV/c)", "ecore >7.0 (GeV/c)"
  };

  TFile *file = new TFile("/direct/phenix+hhj/arbint/RootFiles/CuAuEMCal.root", "r");

  //Get the 3D histogram
  TH3D *hEcore[NARMSECT];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      char hname[256];
      sprintf(hname, "hEmc3D_Armsect_%i", ias);
      hEcore[ias] = (TH3D*)file->Get(hname);
    }


  unsigned int nmax[NARMSECT][NECOREBIN];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          nmax[ias][iec] = 0;
          for (int iz = 0; iz < ZPOS_PBGL; iz++ )
            {
              for (int iy = 0; iy < YPOS_PBGL; iy++ )
                {
                  if (IsValidYZ(ias, iy, iz))
                    {
                      unsigned int nhit = 0;
                      for (int iecc = iec; iecc < NECOREBIN; iecc++)
                        {
                          nhit += hEcore[ias]->GetBinContent(hEcore[ias]->GetXaxis()->FindBin(iz),
                                                             hEcore[ias]->GetYaxis()->FindBin(iy),
                                                             iecc + 1);
			}
                      if (nhit > nmax[ias][iec]) nmax[ias][iec] = nhit;
                    }
                }
            }
        }
    }

  //Make 1D histograms
  TH1D *hHit[NARMSECT][NECOREBIN];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          char hname[256];
          sprintf(hname, "hHit_Armsect%i_ecore%i", ias, iec);

          char htitle[256];
          sprintf(htitle, "%s, %s", armsect[ias], ecore[iec]);

	  int nbins = -1;
          int iass = -1;
	  if (ias == 4 || ias == 5){
	    iass = 1;
	    nbins = 17;
	  }else{
	    iass = 0;
	    nbins = 36;
	  }

	  hHit[ias][iec] = new TH1D(hname, "", nbins, 0, nmax[iass][iec]*0.75);
	  hHit[ias][iec]->SetMarkerStyle(24);
          hHit[ias][iec]->SetMarkerColor(1);
          hHit[ias][iec]->SetMarkerSize(0.6);
          hHit[ias][iec]->SetTitle(htitle);
          hHit[ias][iec]->GetXaxis()->SetTitle("Total Hits");
          hHit[ias][iec]->GetYaxis()->SetTitle("Number of Towers");
	  hHit[ias][iec]->GetYaxis()->SetTitleOffset(1.5);
          TGaxis::SetMaxDigits(4);

          for (int iz = 0; iz < ZPOS_PBGL; iz++ )
            {
              for (int iy = 0; iy < YPOS_PBGL; iy++ )
                {
                  if (IsValidYZ(ias, iy, iz))
                    {
                      float nhitt = 0.0;
                      for (int iecc = iec; iecc < NECOREBIN; iecc++)
                        {
                          nhitt += hEcore[ias]->GetBinContent(hEcore[ias]->GetXaxis()->FindBin(iz),
                                                             hEcore[ias]->GetYaxis()->FindBin(iy),
                                                             iecc + 1);
                        }
                      hHit[ias][iec]->Fill(nhitt);
		      }
                    }
                }
            }
    }


  //Fit the histograms
  float mean[NARMSECT][NECOREBIN];
  float sigma[NARMSECT][NECOREBIN];

  float p0[NARMSECT][NECOREBIN];
  float p1[NARMSECT][NECOREBIN];
  float p2[NARMSECT][NECOREBIN];
  float p3[NARMSECT][NECOREBIN];
  float p4[NARMSECT][NECOREBIN];
  float p5[NARMSECT][NECOREBIN];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
	  hHit[ias][iec]->Sumw2();

          int iass = 0;
          if (ias == 4 || ias == 5) iass = 1;

          mean[ias][iec] = 0.0;
          sigma[ias][iec] = 0.0;

	  float meann = hHit[ias][iec]->GetBinCenter(hHit[ias][iec]->GetMaximumBin());
          float rms = hHit[ias][iec]->GetRMS();
          float inv = hHit[ias][iec]->Integral();

          TF1 *fpeak = new TF1("fpeak", "gaus(0)+gaus(3)", 0, nmax[iass][iec]);
          fpeak->SetParameters(inv / 50, meann, rms, inv / 100, meann, rms*8);

	  
	  fpeak->SetParLimits(1, meann - nmax[iass][iec] / 40, meann + nmax[iass][iec] / 40);
          fpeak->SetParLimits(2, rms / 2, rms*5);
          fpeak->SetParLimits(4, meann - nmax[iass][iec] / 10, meann + nmax[iass][iec] / 10);

          if (iec == 0)
            {
	      fpeak->SetParLimits(5, rms, rms*7);
            }
          else
            {
              fpeak->SetParLimits(5, rms*2, rms*7);
            }

          fpeak->SetLineColor(2);
          fpeak->SetLineWidth(2);
          fpeak->SetLineStyle(2);

          hHit[ias][iec]->Fit("fpeak", "MRQ", "", 0, nmax[iass][iec]);
          mean[ias][iec] = fpeak->GetParameter(1);
          sigma[ias][iec] = fpeak->GetParameter(2);

          p0[ias][iec] = fpeak->GetParameter(0);
          p1[ias][iec] = fpeak->GetParameter(1);
          p2[ias][iec] = fpeak->GetParameter(2);
          p3[ias][iec] = fpeak->GetParameter(3);
          p4[ias][iec] = fpeak->GetParameter(4);
          p5[ias][iec] = fpeak->GetParameter(5);
        }
    }


  //Plot the histograms with mean and sigma
  TCanvas *cmean[NARMSECT];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      //      if(!(ias==4 || ias==5)) continue;
      cmean[ias] = new TCanvas(Form("%s", armsect[ias]), Form("%s", armsect[ias]), 1250, 600);
      cmean[ias]->SetBorderMode(0);
      cmean[ias]->SetFillColor(10);
      cmean[ias]->SetFrameFillColor(10);
      cmean[ias]->Divide(4, 2);

      for (int iec = 0; iec < NECOREBIN; iec++)
        {
          cmean[ias]->cd(iec + 1);

          TPaveText *t0 = new TPaveText(0.58, 0.75, 0.88, 0.88, "brNDC");
          t0->AddText(Form("#splitline{Mean: %.2e}{Sigma: %.2e}", mean[ias][iec], sigma[ias][iec]));
          t0->SetBorderSize(0);
          t0->SetTextFont(42);
          t0->SetTextSize(0.04);
          t0->SetTextColor(2);
          t0->SetFillColor(10);

          int iass = 0;
          if (ias == 4 || ias == 5) iass = 1;

          //Set signal Gaussian parameters
          TF1 *fsg = new TF1("fsg", "gaus(0)", 0, nmax[iass][iec]);
          fsg->SetParameters(p0[ias][iec], p1[ias][iec], p2[ias][iec]);
          fsg->SetLineWidth(2);
          fsg->SetLineStyle(2);
          fsg->SetLineColor(7);

          //Set background Gaussian parameters
          TF1 *fbg = new TF1("fbg", "gaus(0)", 0, nmax[iass][iec]);
          fbg->SetParameters(p3[ias][iec], p4[ias][iec], p5[ias][iec]);
          fbg->SetLineWidth(2);
          fbg->SetLineStyle(2);
          fbg->SetLineColor(3);

          hHit[ias][iec]->Draw("p");
          fsg->Draw("same");
          fbg->Draw("same");
          t0->Draw();

          cmean[ias]->Update();
        }
      cmean[ias]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/JetAnalyzer/WarnMaps/gaussian_%d.png", ias));
    }


  cout << "Writing mean and sigma" << endl;

  ofstream myfileMean ("meanSigma.txt");
  if (myfileMean.is_open())
    {
      for (int ias = 0; ias < NARMSECT; ias++)
        {
          for (int iec = 0; iec < NECOREBIN; iec++)
            {
              myfileMean << mean[ias][iec] << " "
                         << sigma[ias][iec] << "\n";
            }
        }
    }
  myfileMean.close();
}

