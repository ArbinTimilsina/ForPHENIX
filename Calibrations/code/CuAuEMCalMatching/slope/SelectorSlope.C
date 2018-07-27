#define SelectorSlope_cxx

#include "SelectorSlope.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

double cutlow[NPT] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0};
double cuthigh[NPT] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};

int SelectorSlope::GetPtBin(float pT)
{
  int kbin = -1;
  for (int m1 = 0; m1 < NPT; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

void SelectorSlope::Begin(TTree * /*tree*/)
{
  cout << "Begin called" << endl << endl;

  TString option = GetOption();
  rootfname = "/direct/phenix+hhj/arbint/Calibration/rootfiles/slope.root";

  eventcounter = 0;

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hEmcBefore[isect][ipt] = new TH2F(Form("hEmcBefore%d_%d", isect, ipt), Form("Before correction. Sector: %d, p_{T} bin: %d", isect, ipt), 80, -0.4, 0.4, 80, -8.0, 8.0);
          hEmcBefore[isect][ipt]->SetStats(0);
          hEmcBefore[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hEmcBefore[isect][ipt]->GetYaxis()->SetTitle("emcdz (cm)");

          hEmcAfterMean[isect][ipt] = new TH2F(Form("hEmcAfterMean%d_%d", isect, ipt), Form("After mean correction. Sector: %d, p_{T} bin: %d", isect, ipt), 80, -0.4, 0.4, 80, -8.0, 8.0);
          hEmcAfterMean[isect][ipt]->SetStats(0);
          hEmcAfterMean[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hEmcAfterMean[isect][ipt]->GetYaxis()->SetTitle("emcdz (cm)");

          hEmcCorrected[isect][ipt] = new TH2F(Form("hEmcCorrected%d_%d", isect, ipt), Form("After correction. Sector: %d, p_{T} bin: %d", isect, ipt), 80, -0.4, 0.4, 80, -8.0, 8.0);
          hEmcCorrected[isect][ipt]->SetStats(0);
          hEmcCorrected[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hEmcCorrected[isect][ipt]->GetYaxis()->SetTitle("emcdz (cm)");

          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              hEmcdZ[isect][ipt][itheta] = new TH1F(Form("hEmcdZ%d_%d_%d", isect, ipt, itheta), Form("hEmcdZ%d_%d_%d", isect, ipt, itheta), 200, -25.00, 25.00);
            }
        }
    }

  //Load the mean
  LoadMean("mean.txt");

  //Load the sigma
  LoadSigma("sigma.txt");
}

void SelectorSlope::SlaveBegin(TTree * /*tree*/)
{
  cout << "SlaveBegin called" << endl << endl;

  TString option = GetOption();

}

Bool_t SelectorSlope::Process(Long64_t entry)
{
  b_armsect->GetEntry(entry);
  b_zed->GetEntry(entry);
  b_centrality->GetEntry(entry);
  b_pT->GetEntry(entry);
  b_charge->GetEntry(entry);
  b_phi->GetEntry(entry);
  b_emcdphi->GetEntry(entry);
  b_emcdz->GetEntry(entry);
  b_theta->GetEntry(entry);
  b_beta->GetEntry(entry);

  if (eventcounter % 100000000 == 0)
    {
      cout << "Processed " << eventcounter << " events" << endl;
    }

  eventcounter++;

  int armsects = (int)armsect;

  int ptBin = GetPtBin(pT);

  float tanTheta = tan((3.1416 / 2) - beta);

  float mean = mean0[armsects][ptBin] + mean1[armsects][ptBin]*tanTheta +  mean2[armsects][ptBin]*tanTheta*tanTheta*tanTheta;
  float sigma = sigma0[armsects][ptBin] + sigma0[armsects][ptBin]*tanTheta*tanTheta + sigma1[armsects][ptBin]*tanTheta*tanTheta*tanTheta*tanTheta;

  float emcdz_corrected = (emcdz - mean)/(sigma);

  hEmcBefore[armsects][ptBin]->Fill(tanTheta, emcdz);
  hEmcAfterMean[armsects][ptBin]->Fill(tanTheta, emcdz - mean);
  hEmcCorrected[armsects][ptBin]->Fill(tanTheta, emcdz_corrected);

  int thetaBin;
  if (tanTheta > -0.4 && tanTheta <= -0.3) thetaBin = 0;
  if (tanTheta > -0.3 && tanTheta <= -0.2) thetaBin = 1;
  if (tanTheta > -0.2 && tanTheta <= -0.1) thetaBin = 2;
  if (tanTheta > -0.1 && tanTheta <= 0.0) thetaBin = 3;
  if (tanTheta > 0.0 && tanTheta <= 0.1) thetaBin = 4;
  if (tanTheta > 0.1 && tanTheta <= 0.2) thetaBin = 5;
  if (tanTheta > 0.2 && tanTheta <= 0.3) thetaBin = 6;
  if (tanTheta > 0.3 && tanTheta <= 0.4) thetaBin = 7;

  hEmcdZ[armsects][ptBin][thetaBin]->Fill(emcdz);

  //Do the offset correction


  return kTRUE;
}

void SelectorSlope::SlaveTerminate()
{
  cout << "SlaveTerminate called" << endl << endl;

}

void SelectorSlope::Terminate()
{
  cout << "Terminate called" << endl << endl;

  TFile *out = new TFile(rootfname.c_str(), "RECREATE");

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hEmcBefore[isect][ipt]->Write();
          hEmcAfterMean[isect][ipt]->Write();
          hEmcCorrected[isect][ipt]->Write();
          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              hEmcdZ[isect][ipt][itheta]->Write();
            }
        }
    }
  out->Close();
}

void SelectorSlope::LoadMean(const char* filename)
{
  cout << "Loading mean from file " << filename << endl;

  ifstream fmean;
  fmean.open(filename);

  if (!fmean)
    {
      cout << "ERROR!!! File " << filename << " could not be located!" << endl;
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
  cout << "Successfully loaded mean from  file " << filename << endl << endl;

}

void SelectorSlope::LoadSigma(const char* filename)
{
  cout << "Loading sigma from file " << filename << endl;

  ifstream fsigma;
  fsigma.open(filename);

  if (!fsigma)
    {
      cout << "ERROR!!! File " << filename << " could not be located!" << endl;
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
  cout << "Successfully loaded sigma from  file " << filename << endl << endl;

}
