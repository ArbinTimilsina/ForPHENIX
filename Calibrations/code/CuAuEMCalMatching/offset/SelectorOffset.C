#define SelectorOffset_cxx

#include "SelectorOffset.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

double cutlow[NPT] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0};
double cuthigh[NPT] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 6.0, 10.0, 25.0};

int SelectorOffset::GetPtBin(float pT)
{
  int kbin = -1;
  for (int m1 = 0; m1 < NPT; m1++)
    {
      if (pT >= cutlow[m1] && pT < cuthigh[m1])kbin = m1;
    }
  return kbin;
}

void SelectorOffset::Begin(TTree * /*tree*/)
{
  cout << "Begin called" << endl << endl;

  TString option = GetOption();
  rootfname = "/direct/phenix+hhj/arbint/Calibration/rootfiles/offset.root";

  eventcounter = 0;

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hEmcBefore[isect][ipt] = new TH2F(Form("hEmcBefore%d_%d", isect, ipt), Form("Before correction. Sector: %d, p_{T} bin: %d", isect, ipt), 80, -0.4, 0.4, 80, -8.0, 8.0);
          hEmcBefore[isect][ipt]->SetStats(0);
          hEmcBefore[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hEmcBefore[isect][ipt]->GetYaxis()->SetTitle("emcdz (cm)");

          hEmcAfter[isect][ipt] = new TH2F(Form("hEmcAfter%d_%d", isect, ipt), Form("After correction. Sector: %d, p_{T} bin: %d", isect, ipt), 80, -0.4, 0.4, 80, -8.0, 8.0);
          hEmcAfter[isect][ipt]->SetStats(0);
          hEmcAfter[isect][ipt]->GetXaxis()->SetTitle("tan(#theta)");
          hEmcAfter[isect][ipt]->GetYaxis()->SetTitle("emcdz (cm)");

          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              hEmcdZ[isect][ipt][itheta] = new TH1F(Form("hEmcdZ%d_%d_%d", isect, ipt, itheta), Form("hEmcdZ%d_%d_%d", isect, ipt, itheta), 200, -25.00, 25.00);
            }
        }
    }

  //Load the slope
  LoadSlope("slope.txt");

}


void SelectorOffset::SlaveBegin(TTree * /*tree*/)
{
  cout << "SlaveBegin called" << endl << endl;

  TString option = GetOption();

}

Bool_t SelectorOffset::Process(Long64_t entry)
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

  hEmcBefore[armsects][ptBin]->Fill(tanTheta, emcdz);


  float offset = slope[armsects][ptBin] * tanTheta;
  hEmcAfter[armsects][ptBin]->Fill(tanTheta, emcdz - offset);

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

void SelectorOffset::SlaveTerminate()
{
  cout << "SlaveTerminate called" << endl << endl;

}

void SelectorOffset::Terminate()
{
  cout << "Terminate called" << endl << endl;

  TFile *out = new TFile(rootfname.c_str(), "RECREATE");

  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          hEmcBefore[isect][ipt]->Write();
          hEmcAfter[isect][ipt]->Write();
          for (int itheta = 0; itheta < NTHETA; itheta++)
            {
              hEmcdZ[isect][ipt][itheta]->Write();
            }
        }
    }
  out->Close();
}


void SelectorOffset::LoadSlope(const char* filename)
{
  cout << "Loading slope from file " << filename << endl;

  ifstream fslope;
  fslope.open(filename);

  if (!fslope)
    {
      cout << "ERROR!!! File " << filename << " could not be located!" << endl;
    }
  for (int isect = 0; isect < NSECT; isect++)
    {
      for (int ipt = 0; ipt < NPT; ipt++)
        {
          fslope >> slope[isect][ipt];
        }
    }

  fslope.close();
  cout << "Successfully loaded slope from  file " << filename << endl << endl;

}
