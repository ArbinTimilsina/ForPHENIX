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
#include <math.h>

using namespace std;

const int NCHARGE = 2;
const int NARM = 2;
const int NZED = 10;
const int NCENT = 6;

//pT from 0.4 to 25 ==> 13 bins
const int NPT = 14;

double cutlow[NPT] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0};
double cuthigh[NPT] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};
float mombins[NPT+1] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};

int GetPtBin(float pT)
{
  int kbin = -1;
  for (int m1 = 0; m1 < NPT; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

void fitsdPhi()
{
  gStyle->SetOptStat(0);

  TFile *f = new TFile("/direct/phenix+hhj/arbint/RootFiles/Pc3EmcMatching.root", "r");

  float min = -5.00;
  float max =  5.00;

  const char *charge[2] = {"Negative", "Positive"};

  const char *arm[2] = {"East", "West"};

  const char *centrality[6] = {"Centrality: 0-10%", "Centrality: 10-20%", "Centrality: 20-30%", "Centrality: 30-40%", "Centrality: 40-60%",
                               "Centrality: 60-100%"
  };

  const char *pT[13] = {"0.4 <p_{T} <0.6", "0.6< p_{T} <0.8", "0.8 < p_{T} <1.0", "1.0< p_{T} <1.2", "1.2< p_{T} <1.4",
                        "1.4< p_{T} <1.6", "1.6< p_{T} <1.8", "1.8< p_{T} <2.0", "2.0< p_{T} <2.5", "2.5< p_{T} <3.0", "3.0< p_{T} <6.0",
                        "6.0< p_{T} <10.0", "10.0< p_{T} <25.0"
  };

  const char *zed[10] = {"-75cm< zed <-60cm", "-60cm< zed <-45cm", "-45cm< zed <-30cm", "-30cm< zed <-15cm", "-15cm< zed <0cm", "0cm< zed <15cm",
                         "15cm< zed <30cm", "30cm< zed <45cm", "45cm< zed <60cm", "60cm< zed <75cm"
  };


  //Get the 2F histograms
  TH2F *h2F[NARM][NZED][NCENT];
  for (int iarm = 0; iarm < NARM; iarm++)
    {
      for (int ized = 0; ized < NZED; ized++)
        {
          for (int icent = 0; icent < NCENT; icent++)
            {
              h2F[iarm][ized][icent] = (TH2F*)f->Get(Form("hPc3sdPhiFinal%d_%d_%d", iarm, ized, icent));
	      h2F[iarm][ized][icent]->RebinX(3);
            }
        }
    }

  //Initialize 1F histograms
  TH1F *h1F[NCHARGE][NARM][NZED][NCENT][NPT];

  int nbinx = h2F[0][0][0]->GetNbinsX(); //bins in x axis
  float binmin = h2F[0][0][0]->GetXaxis()->GetXmin();
  float binmax = h2F[0][0][0]->GetXaxis()->GetXmax();

  //Initialize the sdPhi and esdPhi also
  float sdPhi[NCHARGE][NARM][NZED][NCENT][NPT][nbinx];
  float esdPhi[NCHARGE][NARM][NZED][NCENT][NPT][nbinx];

  //In addition, initialize the mean and sigma historams
  TH1F *hMeansdPhi[NCHARGE][NARM][NZED][NCENT];
  TH1F *hSigmasdPhi[NCHARGE][NARM][NZED][NCENT];

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int iarm = 0; iarm < NARM; iarm++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  int col = ized + 1;
                  if (ized == 9)col = 22;

                  hMeansdPhi[icharge][iarm][ized][icent] = new TH1F(Form("hMeansdPhi_%d_%d_%d_%d", icharge, iarm, ized, icent),
								     Form("hMeansdPhi_%d_%d_%d_%d", icharge, iarm, ized, icent),
								     NPT, mombins);


                  hMeansdPhi[icharge][iarm][ized][icent]->SetTitle(Form("PC3 <sdPhi> (final): %s, %s", charge[icharge], centrality[icent]));
                  hMeansdPhi[icharge][iarm][ized][icent]->SetTitleSize(0.6);
                  hMeansdPhi[icharge][iarm][ized][icent]->SetMarkerStyle(28);
                  hMeansdPhi[icharge][iarm][ized][icent]->SetMarkerColor(col);


                  hSigmasdPhi[icharge][iarm][ized][icent] = new TH1F(Form("hSigmasdPhi_%d_%d_%d_%d", icharge, iarm, ized, icent),
								      Form("hSigmasdPhi_%d_%d_%d_%d", icharge, iarm, ized, icent),
								      NPT, mombins);


                  hSigmasdPhi[icharge][iarm][ized][icent]->SetTitle(Form("PC3 #sigma_{sdPhi} (final): %s, %s", charge[icharge], centrality[icent]));
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetTitleSize(0.6);
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetMarkerStyle(29);
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetMarkerColor(col);

                  for (int ipt = 0; ipt < NPT; ipt++) //pt bin
                    {
                      h1F[icharge][iarm][ized][icent][ipt] = new TH1F(Form("h_PC3sdPhi_%d_%d_%d_%d_%d", icharge, iarm, ized, icent, ipt),
								       Form("h_PC3sdPhi_%d_%d_%d_%d_%d", icharge, iarm, ized, icent, ipt),
								       nbinx, binmin, binmax);
                      h1F[icharge][iarm][ized][icent][ipt]->Sumw2();
                      h1F[icharge][iarm][ized][icent][ipt]->SetMarkerStyle(24);
                      h1F[icharge][iarm][ized][icent][ipt]->SetMarkerColor(1);
                      h1F[icharge][iarm][ized][icent][ipt]->SetTitle("PC3sdPhi");
                      h1F[icharge][iarm][ized][icent][ipt]->GetXaxis()->SetTitle("sdPhi");
                      h1F[icharge][iarm][ized][icent][ipt]->GetXaxis()->SetRangeUser(min, max);

                      for (int ix = 0; ix < nbinx; ix++)  // sdPhi bin
                        {
                          sdPhi[icharge][iarm][ized][icent][ipt][ix] = 0.0;
                          esdPhi[icharge][iarm][ized][icent][ipt][ix] = 0.0;
                        }
                    }
                }
            }
        }
    }

  //Get the sdPhi and esdPhi
  for (int iarm = 0; iarm < NARM; iarm++)
    {
      for (int ized = 0; ized < NZED; ized++)
        {
          for (int icent = 0; icent < NCENT; icent++)
            {
              for (int iy = 1; iy <= h2F[iarm][ized][icent]->GetNbinsY(); iy++)  //pt
                {

                  //Get charge
                  int separation = h2F[iarm][ized][icent]->GetNbinsY() / 2 + 1;
                  int charge = -1;
                  if (iy < separation )
                    {
                      charge = 0;
                    }
                  else
                    {
                      charge = 1;
                    }

                  float pTT = fabs(h2F[iarm][ized][icent]->GetYaxis()->GetBinCenter(iy));
                  int kbin = GetPtBin(pTT);
                  if (kbin < 0)continue;

                  for (int ix = 1; ix <= nbinx; ix++)  //sdPhi
                    {
		      sdPhi [charge][iarm][ized][icent][kbin][ix-1] += h2F[iarm][ized][icent]->GetBinContent(ix, iy);
                      esdPhi[charge][iarm][ized][icent][kbin][ix-1] += pow(h2F[iarm][ized][icent]->GetBinError(ix, iy), 2);
                    }
                }
            }
        }
    }

  //Fill the 1F histogram
  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int iarm = 0; iarm < NARM; iarm++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  for (int ipt = 0; ipt < NPT; ipt++) //pt
                    {
                      for (int ix = 0; ix < nbinx; ix++)  // sdPhi bin (cout some statements here)
                        {
                          h1F[icharge][iarm][ized][icent][ipt]->SetBinContent(ix + 1, sdPhi[icharge][iarm][ized][icent][ipt][ix]);
                          h1F[icharge][iarm][ized][icent][ipt]->SetBinError(ix + 1, sqrt(esdPhi[icharge][iarm][ized][icent][ipt][ix]));
                        }
                    }
                }
            }
        }
    }

  //Fit the histogram
  float mmin = -2.0;
  float mmax = 2.0;

  float p0[NCHARGE][NARM][NZED][NCENT][NPT];
  float p1[NCHARGE][NARM][NZED][NCENT][NPT];
  float p2[NCHARGE][NARM][NZED][NCENT][NPT];
  float p3[NCHARGE][NARM][NZED][NCENT][NPT];
  float p4[NCHARGE][NARM][NZED][NCENT][NPT];
  float p5[NCHARGE][NARM][NZED][NCENT][NPT];
  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int iarm = 0; iarm < NARM; iarm++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  for (int ipt = 1; ipt < NPT; ipt++) //pt- changed here for low bin
                    {
                      float inv = h1F[icharge][iarm][ized][icent][ipt]->Integral();

		      TF1 *fpeak;
		      float mean  = h1F[icharge][iarm][ized][icent][ipt]->GetBinCenter(h1F[icharge][iarm][ized][icent][ipt]->GetMaximumBin());
                      float sigma = 1.0;

                      if(ipt > 10){
                        fpeak = new TF1("fpeak", "gaus(0) + gaus(3)", min, max);
			fpeak->SetParameters(inv / 10, mean, sigma, inv / 15, mean, 4.0);
			h1F[icharge][iarm][ized][icent][ipt]->Rebin(3);
                      }else{
                        fpeak = new TF1("fpeak", "gaus(0) + gaus(3)", min, max);
                        fpeak->SetParameters(inv / 10, mean, sigma, inv / 15, mean, 4.0);
                      }

		      fpeak->SetLineColor(2);
		      fpeak->SetLineWidth(2);
                      fpeak->SetLineStyle(2);

                      fpeak->SetParLimits(1, mean - 0.2, mean + 0.2);
                      fpeak->SetParLimits(2, 0.8, 1.2);
                      fpeak->SetParLimits(4, mean + mmin, mean + mmax);
                      fpeak->SetParLimits(5, 3.0, 6.0);

                      h1F[icharge][iarm][ized][icent][ipt]->Fit("fpeak", "MRQ", "", min, max);

                      //Fill the mean and sigma historams
                      hMeansdPhi[icharge][iarm][ized][icent]->SetBinContent (ipt + 1, fpeak->GetParameter(1));
                      hMeansdPhi[icharge][iarm][ized][icent]->SetBinError (ipt + 1, fpeak->GetParError(1));

                      hSigmasdPhi[icharge][iarm][ized][icent]->SetBinContent (ipt + 1, fpeak->GetParameter(2));
                      hSigmasdPhi[icharge][iarm][ized][icent]->SetBinError (ipt + 1, fpeak->GetParError(2));

		      p0[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(0);
                      p1[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(1);
                      p2[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(2);
                      p3[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(3);
                      p4[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(4);
                      p5[icharge][iarm][ized][icent][ipt] = fpeak->GetParameter(5);
                    }
                }
            }
        }
    }


  //Plot a sample distribution of 1F histograms and fits
  TCanvas *cfit[NCHARGE][NARM][NZED][NCENT];

  for (int icharge = 0; icharge < NCHARGE; icharge++)
    {
      for (int iarm = 0; iarm < NARM; iarm++)
        {
          for (int ized = 0; ized < NZED; ized++)
            {
              for (int icent = 0; icent < NCENT; icent++)
                {
                  //Make a choice here- Drawing all is insane
		  if (!(iarm == 0 && ized == 5 && icharge==0 && icent==0)) continue;
                  cfit[icharge][iarm][ized][icent] = new TCanvas(Form("cfit_%s_%s_%s_%s", charge[icharge], arm[iarm], zed[ized], centrality[icent]),
								  Form("cfit_%s_%s_%s_%s", charge[icharge], arm[iarm], zed[ized], centrality[icent]), 1250, 350);
                  cfit[icharge][iarm][ized][icent]->Divide(7, 2);
                  for (int ipt = 0; ipt < NPT; ipt++) //pt
                    {
                      //Cosmetics to the 1F hisgogram
                      TPaveText *t0 = new TPaveText(0.13, 0.68, 0.43, 0.81, "brNDC");
                      t0->AddText(Form("#splitline{%s Charge}{%s}", charge[icharge], arm[iarm]));
                      t0->SetBorderSize(0);
                      t0->SetTextFont(42);
                      t0->SetTextSize(0.04);
                      t0->SetTextColor(2);
                      t0->SetFillColor(10);

                      TPaveText *t00 = new TPaveText(0.18, 0.53, 0.42, 0.64, "brNDC");
                      t00->AddText(Form("#splitline{%s}{%s}", pT[ipt], ""));
                      t00->SetBorderSize(0);
                      t00->SetTextFont(42);
                      t00->SetTextSize(0.05);
                      t00->SetTextColor(2);
                      t00->SetFillColor(10);

                      TPaveText *t1 = new TPaveText(0.58, 0.68, 0.88, 0.81, "brNDC");
                      t1->AddText(Form("#splitline{%s}{%s}", zed[ized], centrality[icent]));
                      t1->SetBorderSize(0);
                      t1->SetTextFont(42);
                      t1->SetTextSize(0.04);
                      t1->SetTextColor(1);
                      t1->SetFillColor(10);

                      //Set signal Gaussian parameters
                      TF1 *fsg = new TF1("fsg", "gaus(0)", min, max);
                      fsg->SetParameters(p0[icharge][iarm][ized][icent][ipt], p1[icharge][iarm][ized][icent][ipt], p2[icharge][iarm][ized][icent][ipt]);
                      fsg->SetLineWidth(2);
                      fsg->SetLineStyle(2);
                      fsg->SetLineColor(7);

                      //Set background Gaussian parameters
                      TF1 *fbg = new TF1("fbg", "gaus(0)", min, max);
                      fbg->SetParameters(p3[icharge][iarm][ized][icent][ipt], p4[icharge][iarm][ized][icent][ipt],p5[icharge][iarm][ized][icent][ipt]);
                      fbg->SetLineWidth(2);
                      fbg->SetLineStyle(2);
                      fbg->SetLineColor(3);


                      //Draw just a single histo
                      if (icharge == 0 && iarm == 0 && ized == 0 && icent == 0 && ipt == 3)
                        {
                          TCanvas *cSinglesdPhi = new TCanvas("cSinglesdPhi", "Single sdPhi", 1000, 700);
                          h1F[icharge][iarm][ized][icent][ipt]->Draw("p");
			  fsg->Draw("same");
                          fbg->Draw("same");
                          t0->Draw();
                          t00->Draw();
                          t1->Draw();
                          cSinglesdPhi->SaveAs("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/ptFit/ptFit_single.png");
                          cSinglesdPhi->Close();
                        }

                      cfit[icharge][iarm][ized][icent]->cd(ipt + 1);
                      h1F[icharge][iarm][ized][icent][ipt]->Draw("p");
		      fsg->Draw("same");
		      fbg->Draw("same");
		      t0->Draw();
		      t00->Draw();
		      t1->Draw();
                      cfit[icharge][iarm][ized][icent]->Update();
                    }
		  cfit[icharge][iarm][ized][icent]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/ptFit/ptFit_Charge%d_Sect%d_Zed%d_Cent%d.png", icharge, iarm, ized, icent));
                  //cfit[icharge][iarm][ized][icent]->Close();
                }
            }
        }
    }

  //Plot the mean distribution
  TCanvas *cmean[NARM];
  TLegend *legend1[NARM];
  TLegend *legend2[NARM];

  for (int iarm = 0; iarm < NARM; iarm++)
    {
      cmean[iarm] = new TCanvas(Form("Mean_%s", arm[iarm]), Form("Mean_%s", arm[iarm]), 1250, 400);
      cmean[iarm]->SetBorderMode(0);
      cmean[iarm]->SetFillColor(10);
      cmean[iarm]->SetFrameFillColor(10);
      cmean[iarm]->Divide(6, 2);

      legend1[iarm] = new TLegend(0.14, 0.69, 0.34, 0.86);
      legend1[iarm]->SetBorderSize(0);
      legend1[iarm]->SetLineColor(1);
      legend1[iarm]->SetLineStyle(1);
      legend1[iarm]->SetLineWidth(1);
      legend1[iarm]->SetFillColor(10);
      legend1[iarm]->SetTextFont(42);
      legend1[iarm]->SetTextSize(0.035);

      legend2[iarm] = new TLegend(0.53, 0.69, 0.77, 0.86);
      legend2[iarm]->SetBorderSize(0);
      legend2[iarm]->SetLineColor(1);
      legend2[iarm]->SetLineStyle(1);
      legend2[iarm]->SetLineWidth(1);
      legend2[iarm]->SetFillColor(10);
      legend2[iarm]->SetTextFont(42);
      legend2[iarm]->SetTextSize(0.035);

      for (int icent = 0; icent < NCENT; icent++)
        {
          for (int icharge = 0; icharge < NCHARGE; icharge++)
            {

              if (icharge == 0)
                {
                  cmean[iarm]->cd(icent + 1);
                  gPad->SetLogx();
                }
              else
                {
                  cmean[iarm]->cd(icent + 7);
                  gPad->SetLogx();
                }
              for (int ized = 0; ized < NZED; ized++)
                {
                  hMeansdPhi[icharge][iarm][ized][icent]->SetStats(0);
                  hMeansdPhi[icharge][iarm][ized][icent]->SetMaximum(5.0);
                  hMeansdPhi[icharge][iarm][ized][icent]->SetMinimum(-5.0);

                  hMeansdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetTitle("<sdPhi>");
                  hMeansdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetLabelSize(0.03);
                  hMeansdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetTitleSize(0.05);

                  hMeansdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                  hMeansdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetLabelSize(0.03);
                  hMeansdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetTitleSize(0.05);


                  if (ized == 0)
                    {
                      hMeansdPhi[icharge][iarm][ized][icent]->Draw();
                    }
                  else
                    {
                      hMeansdPhi[icharge][iarm][ized][icent]->Draw("same");
                    }

                  if (icharge == 0 && icent == 0 && ized < 5) legend1[iarm]->AddEntry(hMeansdPhi[icharge][iarm][ized][icent], Form("%s", zed[ized]), "p");

                  if (icharge == 0 && icent == 0 && ized >= 5) legend2[iarm]->AddEntry(hMeansdPhi[icharge][iarm][ized][icent], Form("%s", zed[ized]), "p");

                  if (icent == 0 && icharge == 0)
                    {
                      legend1[iarm]->Draw("same");
                      legend2[iarm]->Draw("same");
                    }
                }
            }
          cmean[iarm]->Update();
        }
      cmean[iarm]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/meanSigma/mean_Arm%d.png", iarm));
      cmean[iarm]->Close();
    }

  //Draw just a singe histo for Mean
  TCanvas *cSingleMean = new TCanvas("cSingleMean", "Single Mean", 1000, 700);
  for (int ized = 0; ized < NZED; ized++)
    {
      gPad->SetLogx();
      if (ized == 0)
        {
          hMeansdPhi[0][0][ized][0]->Draw();
        }
      else
        {
          hMeansdPhi[0][0][ized][0]->Draw("same");
        }
      legend1[0]->Draw("same");
      legend2[0]->Draw("same");
      cSingleMean->Update();
    }
  cSingleMean->SaveAs("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/meanSigma/mean_Single.png");
  cSingleMean->Close();



  //Plot the sigma distribution
  TCanvas *csigma[NARM];
  TLegend *legend11[NARM];
  TLegend *legend22[NARM];
  for (int iarm = 0; iarm < NARM; iarm++)
    {
      csigma[iarm] = new TCanvas(Form("Sigma_%s", arm[iarm]), Form("Sigma_%s", arm[iarm]), 1250, 400);
      csigma[iarm]->SetBorderMode(0);
      csigma[iarm]->SetFillColor(10);
      csigma[iarm]->SetFrameFillColor(10);
      csigma[iarm]->Divide(6, 2);

      legend11[iarm] = new TLegend(0.14, 0.69, 0.34, 0.86);
      legend11[iarm]->SetBorderSize(0);
      legend11[iarm]->SetLineColor(1);
      legend11[iarm]->SetLineStyle(1);
      legend11[iarm]->SetLineWidth(1);
      legend11[iarm]->SetFillColor(10);
      legend11[iarm]->SetTextFont(42);
      legend11[iarm]->SetTextSize(0.035);

      legend22[iarm] = new TLegend(0.53, 0.69, 0.77, 0.86);
      legend22[iarm]->SetBorderSize(0);
      legend22[iarm]->SetLineColor(1);
      legend22[iarm]->SetLineStyle(1);
      legend22[iarm]->SetLineWidth(1);
      legend22[iarm]->SetFillColor(10);
      legend22[iarm]->SetTextFont(42);
      legend22[iarm]->SetTextSize(0.035);

      for (int icent = 0; icent < NCENT; icent++)
        {
          for (int icharge = 0; icharge < NCHARGE; icharge++)
            {
              if (icharge == 0)
                {
                  csigma[iarm]->cd(icent + 1);
                  gPad->SetLogx();
                }
              else
                {
                  csigma[iarm]->cd(icent + 7);
                  gPad->SetLogx();
                }
              for (int ized = 0; ized < NZED; ized++)
                {
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetStats(0);
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetMaximum(6.0);
                  hSigmasdPhi[icharge][iarm][ized][icent]->SetMinimum(-1.0);

                  hSigmasdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetTitle("#sigma_{sdPhi}");
                  hSigmasdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetLabelSize(0.03);
                  hSigmasdPhi[icharge][iarm][ized][icent]->GetYaxis()->SetTitleSize(0.05);

                  hSigmasdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                  hSigmasdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetLabelSize(0.03);
                  hSigmasdPhi[icharge][iarm][ized][icent]->GetXaxis()->SetTitleSize(0.05);

                  if (ized == 0)
                    {
                      hSigmasdPhi[icharge][iarm][ized][icent]->Draw();
                    }
                  else
                    {
                      hSigmasdPhi[icharge][iarm][ized][icent]->Draw("same");
                    }

                  if (icharge == 0 && icent == 0 && ized < 5) legend11[iarm]->AddEntry(hSigmasdPhi[icharge][iarm][ized][icent], Form("%s", zed[ized]), "p");

                  if (icharge == 0 && icent == 0 && ized >= 5) legend22[iarm]->AddEntry(hSigmasdPhi[icharge][iarm][ized][icent], Form("%s", zed[ized]), "p");

                  if (icent == 0 && icharge == 0)
                    {
                      legend11[iarm]->Draw("same");
                      legend22[iarm]->Draw("same");
                    }
                }
            }
          csigma[iarm]->Update();
        }
      csigma[iarm]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/meanSigma/sigma_Arm%d.png", iarm));
      csigma[iarm]->Close();
    }

  //Draw just a singe histo for Sigma
  TCanvas *cSingleSigma = new TCanvas("cSingleSigma", "Single Sigma", 1000, 700);
  for (int ized = 0; ized < NZED; ized++)
    {
      gPad->SetLogx();
      if (ized == 0)
        {
          hSigmasdPhi[0][0][ized][0]->Draw();
        }
      else
        {
          hSigmasdPhi[0][0][ized][0]->Draw("same");
        }

      legend11[0]->Draw("same");
      legend22[0]->Draw("same");
      cSingleSigma->Update();
    }
  cSingleSigma->SaveAs("/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/dphi/meanSigma/sigma_Single.png");
  cSingleSigma->Close();
}

