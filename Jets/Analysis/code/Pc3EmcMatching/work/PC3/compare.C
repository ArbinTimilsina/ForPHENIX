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
#include <TNtuple.h>
#include <TTree.h>
#include <TCut.h>

using namespace std;

void compare()
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile("/direct/phenix+hhj/arbint/RootFiles/Pc3EmcMatching.root", "r");

  TH1F *pTAll = new TH1F("pTAll", "Comparision", 50, 0, 50);

  TH1F *pTEmc = new TH1F("pTEmc", "Emc matching", 50, 0, 50);

  TH1F *pTPc3 = new TH1F("pTPc3", "Pc3 matching", 50, 0, 50);

  TH1F *pTBoth = new TH1F("pTBoth", "(Emc || Pc3) matching", 50, 0, 50);
 
  TCanvas *c1 = new TCanvas("c1", "Quality board", 1200, 800);
  pTAll = (TH1F*)file->Get("hAll");
  pTEmc = (TH1F*)file->Get("hEmc");
  pTPc3 = (TH1F*)file->Get("hPc3");
  pTBoth = (TH1F*)file->Get("hBoth");
  c1->Clear();

  pTAll->SetStats(0);
  pTAll->SetMinimum(1);
  pTAll->SetLineColor(1);
  //  pTAll->Sumw2();
  //  pTAll->SetMarkerStyle(28);
  pTAll->SetTitle("Comparision");
  pTAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  pTAll->GetYaxis()->SetTitle("Yield");
  pTAll->GetXaxis()->SetRangeUser(0, 24);

  pTEmc->SetStats(0);
  pTEmc->SetMinimum(1);
  pTEmc->SetLineColor(2);
  //pTEmc->Sumw2();
  // pTEmc->SetMarkerStyle(28);

  pTPc3->SetStats(0);
  pTPc3->SetMinimum(1);
  pTPc3->SetLineColor(4);
  //pTPc3->Sumw2();
  //pTPc3->SetMarkerStyle(28);

  pTBoth->SetStats(0);
  pTBoth->SetMinimum(1);
  pTBoth->SetLineColor(7);
  //pTBoth->Sumw2();
  //pTBoth->SetMarkerStyle(28);

  TLegend *legend = new TLegend(0.55, 0.70, 0.87, 0.87);
  gPad->SetLogy();
  legend->SetTextSize(0.03);
  legend->AddEntry(pTAll, "Good Tracks","l");
  legend->AddEntry(pTEmc, "Emc matching","l");
  legend->AddEntry(pTPc3, "Pc3 matching","l");
  legend->AddEntry(pTBoth, "(Emc || Pc3) matching","l");

  pTAll->Draw();
  pTEmc->Draw("SAME");
  pTPc3->Draw("SAME");
  pTBoth->Draw("SAME");
  legend->Draw();

  //  TString fName = "/direct/phenix+hhj/arbint/plots/Pc3EmcMatching/pT.png";
  //c1->Print(fName);
}
