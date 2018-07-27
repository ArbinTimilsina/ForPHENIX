#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void plotClusters()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/He3Jet.root", "READ");;
    TPaveText *species = getSpecies(3);
    species->AddText("Run 14 He3+Au @ 200 GeV, ERT");

    TH1F *hEvents = (TH1F*)file->Get("hEvents");
    float nEvents = hEvents->GetBinContent(2);

    //************************************************************************************
    //************************************************************************************

    TH1F *hClustersTower = (TH1F*)file->Get("hClustersTower");
    hClustersTower->SetMarkerStyle(6);
    hClustersTower->SetMarkerColor(2);

    hClustersTower->Sumw2();
    hClustersTower->Scale(1 / nEvents);
    setHisto(hClustersTower, "Tower id", "Hits per event");
    setHistoMin(hClustersTower, true, 3);
    hClustersTower->SetMaximum(5E-3);
    hClustersTower->SetMinimum(4E-8);

    TPaveText *titleTower = getHistoTitle();
    titleTower->AddText("");

    TCanvas *cClusters2 = new TCanvas("cClusters2", "Clusters going into Jet Reconstruction", 1200, 450);
    cClusters2->Clear();
    cClusters2->cd();
    gPad->SetLogy();
    hClustersTower->Draw();
    titleTower->Draw("SAME");
    species->Draw("SAME");
    cClusters2->Update();
    cClusters2->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/plots/ClustersTower.png");
}




