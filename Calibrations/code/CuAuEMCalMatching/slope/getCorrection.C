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

void getCorrection()
{
  gStyle->SetOptStat(0);

  TFile *f = new TFile("/direct/phenix+hhj/arbint/Calibration/rootfiles/slope.root", "r");

  float min = -13.0;
  float max =  13.0;

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
              h1F[isect][ipt][itheta]->Sumw2();
              h1F[isect][ipt][itheta]->SetMarkerStyle(24);
              h1F[isect][ipt][itheta]->SetMarkerColor(1);
              h1F[isect][ipt][itheta]->SetTitle(Form("hEmcdZ: %s, %s", sector[isect], pT[ipt]));
              h1F[isect][ipt][itheta]->GetXaxis()->SetTitle("EMCdZ");
              h1F[isect][ipt][itheta]->GetXaxis()->SetRangeUser(min, max);

              if (ipt > 10)  h1F[isect][ipt][itheta]->GetXaxis()->SetRangeUser(-6.0, 6.0);
            }
        }
    }

  TH1F *hMeandZ[NSECT][NPT];
  TH1F *hSigmadZ[NSECT][NPT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hMeandZ[isect][ipt] = new TH1F(Form("hMeandZ_%d_%d", isect, ipt),
                                         Form("hMeandZ_%d_%d", isect, ipt),
                                         NTHETA, thetabins);


          hMeandZ[isect][ipt]->SetTitle(Form("EMC mean dz: %s, %s", sector[isect], pT[ipt]));
          hMeandZ[isect][ipt]->Sumw2();
          hMeandZ[isect][ipt]->SetMarkerStyle(28);
          hMeandZ[isect][ipt]->SetMarkerColor(1);
          hMeandZ[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hMeandZ[isect][ipt]->GetYaxis()->SetTitle("<emcdz> cm");
          hMeandZ[isect][ipt]->GetYaxis()->SetRangeUser(-6.0, 6.0);;

          hSigmadZ[isect][ipt] = new TH1F(Form("hSigmadZ_%d_%d", isect, ipt),
                                          Form("hSigmadZ_%d_%d", isect, ipt),
                                          NTHETA, thetabins);


          hSigmadZ[isect][ipt]->SetTitle(Form("EMC sigma dz: %s, %s", sector[isect], pT[ipt]));
          hSigmadZ[isect][ipt]->Sumw2();
          hSigmadZ[isect][ipt]->SetMarkerStyle(28);
          hSigmadZ[isect][ipt]->SetMarkerColor(1);
          hSigmadZ[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hSigmadZ[isect][ipt]->GetYaxis()->SetTitle("#sigma_{emcdz}");
          hSigmadZ[isect][ipt]->GetYaxis()->SetRangeUser(0.0, 8.0);;
        }
    }

  //Fit the histogram
  float p0[NSECT][NPT][NTHETA];
  float p1[NSECT][NPT][NTHETA];
  float p2[NSECT][NPT][NTHETA];
  float p3[NSECT][NPT][NTHETA];
  float p4[NSECT][NPT][NTHETA];
  float p5[NSECT][NPT][NTHETA];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          for (int itheta = 0; itheta < NTHETA; itheta++)
            {

              float inv = h1F[isect][ipt][itheta]->Integral();
              float mean  =  h1F[isect][ipt][itheta]->GetBinCenter(h1F[isect][ipt][itheta]->GetMaximumBin());

              float minn;
              float maxx;

              if (ipt > 10)
                {
                  minn = -6.0;
                  maxx = 6.0;
                }
              else
                {
                  minn = min;
                  maxx = max;
                }

              TF1 *fpeak;
              if (ipt > 10)
                {
                  fpeak = new TF1("fpeak", "gaus(0) + gaus(3)", minn, maxx);
                  fpeak->SetParameters(inv / 10, mean, 1.5, inv / 15, mean, 5.0);
                  fpeak->SetParLimits(1, mean - 2.0, mean + 2.0);
                  fpeak->SetParLimits(2, 1.0, 3.0);
                  fpeak->SetParLimits(4, mean - 2.0, mean + 2.0);
                  fpeak->SetParLimits(5, 4.5, 6.0);
                }
              else if (ipt == 0)
                {
                  fpeak = new TF1("fpeak", "gaus(0) + gaus(3)", minn, maxx);
                  fpeak->SetParameters(inv / 10, mean, 3.0, inv / 15, mean, 15);
                  fpeak->SetParLimits(1, mean - 0.5, mean + 0.5);
                  fpeak->SetParLimits(2, 2.0, 7.0);
                  fpeak->SetParLimits(4, mean - 6.0, mean + 6.0);
                  fpeak->SetParLimits(5, 13, 25.0);
                }
              else
                {
                  fpeak = new TF1("fpeak", "gaus(0) + gaus(3) + pol1(6)", minn, maxx);
                  fpeak->SetParameters(inv / 10, mean, 1.5, inv / 15, mean, 6.0, 0.0, 0.0);
                  fpeak->SetParLimits(1, mean - 0.1, mean + 0.1);
                  fpeak->SetParLimits(2, 1.0, 3.5);
                  fpeak->SetParLimits(4, mean - 6.0, mean + 6.0);
                  fpeak->SetParLimits(5, 4.0, 10.0);
                }

              fpeak->SetLineColor(2);
              fpeak->SetLineWidth(2);
              fpeak->SetLineStyle(2);

              h1F[isect][ipt][itheta]->Fit("fpeak", "MRQ", "", minn, maxx);

              //Fill the mean histogram
              hMeandZ[isect][ipt]->SetBinContent (itheta + 1, fpeak->GetParameter(1));
              hMeandZ[isect][ipt]->SetBinError (itheta + 1, fpeak->GetParError(1));

              //Fill the sigma histogram
              hSigmadZ[isect][ipt]->SetBinContent (itheta + 1, fpeak->GetParameter(2));
              hSigmadZ[isect][ipt]->SetBinError (itheta + 1, fpeak->GetParError(2));

              p0[isect][ipt][itheta] = fpeak->GetParameter(0);
              p1[isect][ipt][itheta] = fpeak->GetParameter(1);
              p2[isect][ipt][itheta] = fpeak->GetParameter(2);
              p3[isect][ipt][itheta] = fpeak->GetParameter(3);
              p4[isect][ipt][itheta] = fpeak->GetParameter(4);
              p5[isect][ipt][itheta] = fpeak->GetParameter(5);
            }
        }
    }

  float mminn = -0.4;
  float mmaxx = 0.4;

  float mean0[NSECT][NPT];
  float mean1[NSECT][NPT];
  float mean2[NSECT][NPT];
  //Fit the mean distribution
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TF1 *fpeak = new TF1("fpeak", "[0]+[1]*x+[2]*x*x*x", mminn, mmaxx);
          fpeak->SetParameters(0.0, 0.0, 0.0);
	  fpeak->SetParLimits(1, -15.0, 15.0);
	  fpeak->SetParLimits(2, -25.0, 25.0);
	 
          fpeak->SetLineColor(2);
          fpeak->SetLineWidth(2);
          fpeak->SetLineStyle(2);

          hMeandZ[isect][ipt]->Fit("fpeak", "MRQ", "", mminn, mmaxx);

          mean0[isect][ipt] = fpeak->GetParameter(0);
          mean1[isect][ipt] = fpeak->GetParameter(1);
          mean2[isect][ipt] = fpeak->GetParameter(2);
        }
    }

  float sigma0[NSECT][NPT];
  float sigma1[NSECT][NPT];
  float sigma2[NSECT][NPT];
  //Fit the sigma distribution
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TF1 *fpeak = new TF1("fpeak", "[0]+[1]*x*x+[2]*x*x*x*x", mminn, mmaxx);
          fpeak->SetParameters(0.0, 0.0, 0.0);
          fpeak->SetParLimits(1, -10.0, 10.0);
          fpeak->SetParLimits(2, -50.0, 50.0);

          fpeak->SetLineColor(2);
          fpeak->SetLineWidth(2);
          fpeak->SetLineStyle(2);

          hSigmadZ[isect][ipt]->Fit("fpeak", "MRQ", "", mminn, mmaxx);

          sigma0[isect][ipt] = fpeak->GetParameter(0);
          sigma1[isect][ipt] = fpeak->GetParameter(1);
          sigma2[isect][ipt] = fpeak->GetParameter(2);
        }
    }

  //Plot the gaussian fit of the histogram
  TCanvas *cfit[NSECT][NPT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
	  //  if (!(isect == 4 && ipt==3)) continue;
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

              float minn;
              float maxx;
              if (ipt > 10)
                {
                  minn = -6.0;
                  maxx = 6.0;
                }
              else
                {
                  minn = min;
                  maxx = max;
                }

              //Set signal Gaussian parameters
              TF1 *fsg = new TF1("fsg", "gaus(0)", minn, maxx);
              fsg->SetParameters(p0[isect][ipt][itheta], p1[isect][ipt][itheta], p2[isect][ipt][itheta]);
              fsg->SetLineWidth(2);
              fsg->SetLineStyle(2);
              fsg->SetLineColor(7);

              //Set background Gaussian parameters
              TF1 *fbg = new TF1("fbg", "gaus(0)", minn, maxx);
              fbg->SetParameters(p3[isect][ipt][itheta], p4[isect][ipt][itheta], p5[isect][ipt][itheta]);
              fbg->SetLineWidth(2);
              fbg->SetLineStyle(2);
              fbg->SetLineColor(3);

              h1F[isect][ipt][itheta]->Draw("p");
              fbg->Draw("lsame");
              fsg->Draw("lsame");

              t0->Draw();
              cfit[isect][ipt]->Update();
            }
          cfit[isect][ipt]->SaveAs(Form("/direct/phenix+hhj/arbint/Calibration/plots/dz/correction/gaussianFit/sector%d_pt%d.png", isect, ipt));
          cfit[isect][ipt]->Close();
        }
    }

  //Plot the mean distribution
  TCanvas *cmean[NSECT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      //      if (!(isect == 0)) continue;
      cmean[isect] = new TCanvas(Form("Mean_%s", sector[isect]), Form("Mean_%s", sector[isect]), 1250, 400);
      cmean[isect]->SetBorderMode(0);
      cmean[isect]->SetFillColor(10);
      cmean[isect]->SetFrameFillColor(10);
      cmean[isect]->Divide(7, 2);

      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TPaveText *t0 = new TPaveText(0.13, 0.68, 0.43, 0.81, "brNDC");
          t0->AddText(Form("%s", pT[ipt]));
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
      cmean[isect]->SaveAs(Form("/direct/phenix+hhj/arbint/Calibration/plots/dz/correction/mean/MeanSector%d.png", isect));
      cmean[isect]->Close();
    }


  //Plot the sigma distribution
  TCanvas *csigma[NSECT];
  for (int isect = 0; isect < NSECT; isect++)
    {
      //      if (!(isect == 2)) continue;
      csigma[isect] = new TCanvas(Form("Sigma_%s", sector[isect]), Form("Sigma_%s", sector[isect]), 1250, 400);
      csigma[isect]->SetBorderMode(0);
      csigma[isect]->SetFillColor(10);
      csigma[isect]->SetFrameFillColor(10);
      csigma[isect]->Divide(7, 2);

      for (int ipt = 0; ipt < NPT; ipt++)
        {
          TPaveText *t0 = new TPaveText(0.13, 0.68, 0.43, 0.81, "brNDC");
          t0->AddText(Form("%s", pT[ipt]));
          t0->SetBorderSize(0);
          t0->SetTextFont(42);
          t0->SetTextSize(0.04);
          t0->SetTextColor(2);
          t0->SetFillColor(10);

          csigma[isect]->cd(ipt + 1);
          hSigmadZ[isect][ipt]->Draw("p");
          t0->Draw("same");
          csigma[isect]->Update();
        }
      csigma[isect]->SaveAs(Form("/direct/phenix+hhj/arbint/Calibration/plots/dz/correction/sigma/SigmaSector%d.png", isect));
      csigma[isect]->Close();
    }


  //Write out a text file
  cout << "Writing out the mean in a text file" << endl;

  ofstream myfileMean("mean.txt");
  if (myfileMean.is_open())
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ipt = 0; ipt < NPT; ipt++)
            {
              myfileMean << mean0[isect][ipt] << " "
			 << mean1[isect][ipt] << " "
			 << mean2[isect][ipt] << "\n";
            }
        }
    }
  myfileMean.close();

  cout << "Writing out the sigma in a text file" << endl;
  ofstream myfileSigma ("sigma.txt");
  if (myfileSigma.is_open())
    {
      for (int isect = 0; isect < NSECT; isect++)
        {
          for (int ipt = 0; ipt < NPT; ipt++)
            {
              myfileSigma << sigma0[isect][ipt] << " "
              << sigma1[isect][ipt] << " "
              << sigma2[isect][ipt] << "\n";
            }
        }
    }
  myfileSigma.close();
}
