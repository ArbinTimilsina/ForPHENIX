#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>

using namespace std;

const int NSECT = 8;
const int NPT = 14;
const int NTHETA = 8;

float mombins[NPT+1] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};
float thetabins[NTHETA+1] = { -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4};

void getSlope()
{
  gStyle->SetOptStat(0);

  TFile *f = new TFile("/direct/phenix+hhj/arbint/Calibration/rootfiles/offset.root", "r");

  float min = -25.0;
  float max =  25.0;

  const char *sector[8] = {"Sector: 0(PbSc)", "Sector: 1(PbSc)", "Sector: 2(PbSc)", "Sector: 3(PbSc)", "Sector: 4(PbGl)", "Sector: 5(PbGl)",
                           "Sector: 6(PbSc)", "Sector: 7(PbSc)"
  };

  const char *pT[14] = {"0.2< p_{T} <0.4", "0.4 <p_{T} <0.6", "0.6< p_{T} <0.8", "0.8 < p_{T} <1.0", "1.0< p_{T} <1.2", "1.2< p_{T} <1.4",
                        "1.4< p_{T} <1.6", "1.6< p_{T} <1.8", "1.8< p_{T} <2.0", "2.0< p_{T} <2.5", "2.5< p_{T} <3.0", "3.0< p_{T} <6.0",
                        "6.0< p_{T} <10.0", "10.0< p_{T} <25.0"
  };

  const char *tTheta[8] = {"-0.4< tan(#theta) <-0.3", "-0.3< tan(#theta) <-0.2", "-0.2< tan(#theta) <-0.1", "-0.1< tan(#theta) <0.0",
                           "0.0< tan(#theta) <0.1", "0.1< tan(#theta) <0.2", "0.2< tan(#theta) <0.3", "0.3< tan(#theta) <0.4"

  };

  TH1F *h1F[NSECT][NPT][NTHETA];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              h1F[isect][ipt][itheta] = (TH1F*)f->Get(Form("hEmcdZ%d_%d_%d", isect, ipt, itheta));
              h1F[isect][ipt][itheta]->Rebin(4);
              h1F[isect][ipt][itheta]->Sumw2();
              h1F[isect][ipt][itheta]->SetMarkerStyle(24);
              h1F[isect][ipt][itheta]->SetMarkerColor(1);
              h1F[isect][ipt][itheta]->SetTitle(Form("hEmcdZ: %s, %s", sector[isect], pT[ipt]));
              h1F[isect][ipt][itheta]->GetXaxis()->SetTitle("EMCdZ");

            }
        }
    }

  TH1F *hMeandZ[NSECT][NPT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hMeandZ[isect][ipt] = new TH1F(Form("hMeandZ_%d_%d", isect, ipt),
                                         Form("hMeandZ_%d_%d", isect, ipt),
                                         NTHETA, thetabins);


          hMeandZ[isect][ipt]->SetTitle(Form("EMC dz: %s, %s", sector[isect], pT[ipt]));
          hMeandZ[isect][ipt]->Sumw2();
          hMeandZ[isect][ipt]->SetMarkerStyle(28);
          hMeandZ[isect][ipt]->SetMarkerColor(1);
          hMeandZ[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hMeandZ[isect][ipt]->GetYaxis()->SetTitle("<emcdz> cm");
          hMeandZ[isect][ipt]->GetYaxis()->SetRangeUser(-6.0, 6.0);;
        }
    }

  //Fit the histogram
  float mmin = -5.00;
  float mmax = 5.00;

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          for (int itheta = 0; itheta < NTHETA; itheta++)
            {

              float inv = h1F[isect][ipt][itheta]->Integral();
              float mean  =  h1F[isect][ipt][itheta]->GetBinCenter(h1F[isect][ipt][itheta]->GetMaximumBin());
              float sigma = 2.5;

              TF1 *fpeak;
              if (ipt > 10)
                {
                  fpeak = new TF1("fpeak", "gaus(0) + gaus(3)", min, max);
                  fpeak->SetParameters(inv / 10, mean, sigma, inv / 15, mean, 12.0);
                }
              else
                {
                  fpeak = new TF1("fpeak", "gaus(0) + gaus(3) + pol1(6)", min, max);
                  fpeak->SetParameters(inv / 10, mean, sigma, inv / 15, mean, 12.0, 0.0, 0.0);
                }

              fpeak->SetLineColor(2);
              fpeak->SetLineWidth(2);
              fpeak->SetLineStyle(2);

              fpeak->SetParLimits(1, mean + mmin, mean + mmax);
              fpeak->SetParLimits(2, 1, 5);
              fpeak->SetParLimits(4, mean + mmin, mean + mmax);
              fpeak->SetParLimits(5, 9, 18);

              h1F[isect][ipt][itheta]->Fit("fpeak", "MRQ", "", min, max);

              //Fill the mean histogram
              hMeandZ[isect][ipt]->SetBinContent (itheta + 1, fpeak->GetParameter(1));
              hMeandZ[isect][ipt]->SetBinError (itheta + 1, fpeak->GetParError(1));
            }
        }
    }


  TH1F *hOffset[NSECT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      hOffset[isect] = new TH1F(Form("hOffset_%d", isect),
                                Form("hOffset_%d", isect),
                                NPT, mombins);

      hOffset[isect]->SetTitle(Form("Offset: %s", sector[isect]));
      hOffset[isect]->Sumw2();
      hOffset[isect]->SetMarkerStyle(28);
      hOffset[isect]->SetMarkerColor(1);
      hOffset[isect]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hOffset[isect]->GetYaxis()->SetTitle("Slope");
    }



  float minn = -0.4;
  float maxx = 0.4;
  float slope[NSECT][NPT];

  //Fit the mean distribution
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TF1 *fpeak = new TF1("fpeak", "[0]+x*[1]", minn, maxx);
          fpeak->SetLineColor(2);
          fpeak->SetLineWidth(2);
          fpeak->SetLineStyle(2);

          hMeandZ[isect][ipt]->Fit("fpeak", "MRQ", "", minn, maxx);

          slope[isect][ipt] = fpeak->GetParameter(1);

          //Fill the histogram
          hOffset[isect]->SetBinContent (ipt + 1, fpeak->GetParameter(1));
          hOffset[isect]->SetBinError (ipt + 1, fpeak->GetParError(1));
        }
    }

  //Plot the gaussian fit of the histogram
  TCanvas *cfit[NSECT][NPT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
	  //          if (!(isect == 0 && ipt == 0)) continue;

          cfit[isect][ipt] = new TCanvas(Form("cfit_%d_%d", isect, ipt), Form("cfit_%s_%s", sector[isect], pT[ipt]), 1250, 450);
          cfit[isect][ipt]->Divide(4, 2);

          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              TPaveText *t0 = new TPaveText(0.13, 0.68, 0.43, 0.81, "brNDC");
              t0->AddText(Form("%s", tTheta[itheta]));
              t0->SetBorderSize(0);
              t0->SetTextFont(42);
              t0->SetTextSize(0.04);
              t0->SetTextColor(2);
              t0->SetFillColor(10);

              cfit[isect][ipt]->cd(itheta + 1);
              h1F[isect][ipt][itheta]->Draw("p");
              t0->Draw("same");
              cfit[isect][ipt]->Update();
            }
          cfit[isect][ipt]->SaveAs(Form("/direct/phenix+hhj/arbint/Calibration/plots/offset/gaussianFit/sector%d_pt%d.png", isect, ipt));
          cfit[isect][ipt]->Close();
        }
    }

  //Plot the mean distribution
  TCanvas *cmean[NSECT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      //if (!(isect == 0)) continue;
      cmean[isect] = new TCanvas(Form("Mean_%s", sector[isect]), Form("Mean_%s", sector[isect]), 1250, 400);
      cmean[isect]->SetBorderMode(0);
      cmean[isect]->SetFillColor(10);
      cmean[isect]->SetFrameFillColor(10);
      cmean[isect]->Divide(7, 2);

      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TPaveText *t0 = new TPaveText(0.13, 0.68, 0.43, 0.81, "brNDC");
          t0->AddText(Form("#splitline{%s}{Slope is: %.2f}", pT[ipt], slope[isect][ipt]));
          t0->SetBorderSize(0);
          t0->SetTextFont(42);
          t0->SetTextSize(0.04);
          t0->SetTextColor(2);
          t0->SetFillColor(10);

          cmean[isect]->cd(ipt + 1);
          hMeandZ[isect][ipt]->Draw("p");
          t0->Draw("same");
          cmean[isect]->Update();
        }
      cmean[isect]->SaveAs(Form("/direct/phenix+hhj/arbint/Calibration/plots/offset/mean/sector%d.png", isect));
      cmean[isect]->Close();
    }

  //Write out a text file
  cout<<"Writing out the slopes in a text file"<<endl;

  ofstream myfileSlope ("slope.txt");
  if (myfileSlope.is_open())
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ipt = 0; ipt < NPT; ipt++)
            {
              myfileSlope << slope[isect][ipt] << "\n";
            }
        }
    }
  myfileSlope.close();
}
