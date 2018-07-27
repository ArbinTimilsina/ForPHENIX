#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void plotPc3()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/pAuJet/pAuJet.root", "READ");;
    TPaveText *species = getSpecies(3);
    species->AddText("Run 15 p+Au @ 200 GeV, ERT");

    //***************************************************************************************
    //***************************************************************************************

    TH1F *hPC3dphi = (TH1F*)file->Get("hPC3dphi");
    hPC3dphi->Rebin(6);
    hPC3dphi->GetXaxis()->SetRangeUser(-0.02, 0.02);
    hPC3dphi->SetLineColor(1);
    hPC3dphi->SetLineWidth(2);
    setHisto(hPC3dphi, "pc3dphi", "Counts/Bin");
    setHistoMax(hPC3dphi, false, 0.2);

    float meanPC3dphi = hPC3dphi->GetBinCenter(hPC3dphi->GetMaximumBin());
    float rmsPC3dphi = hPC3dphi->GetRMS();
    float invPC3dphi = hPC3dphi->Integral();

    TF1 *fPeakPC3dphi = new TF1("fPeakPC3dphi", "gaus(0)", -0.0045, 0.005);
    fPeakPC3dphi->SetParameters(invPC3dphi / 50, meanPC3dphi, rmsPC3dphi);
    fPeakPC3dphi->SetLineColor(2);

    hPC3dphi->Fit("fPeakPC3dphi", "MRQ", "", -0.0045, 0.005);

    TPaveText *tPeakPC3dphi = new TPaveText(0.435, 0.774, 0.886, 0.848, "brNDC");
    tPeakPC3dphi->SetBorderSize(0);
    tPeakPC3dphi->SetTextFont(2);
    tPeakPC3dphi->SetTextSize(0.04);
    tPeakPC3dphi->SetTextColor(1);
    tPeakPC3dphi->SetFillColor(10);
    tPeakPC3dphi->AddText(Form("Mean=%.3e, Sigma=%.3e", fPeakPC3dphi->GetParameter(1), fPeakPC3dphi->GetParameter(2)));

    TH1F *hPC3dz = (TH1F*)file->Get("hPC3dz");
    hPC3dz->GetXaxis()->SetRangeUser(-15.0, 15.0);
    hPC3dz->SetLineColor(1);
    hPC3dz->SetLineWidth(2);
    setHisto(hPC3dz, "pc3dz", "Counts/Bin");
    setHistoMax(hPC3dz, false, 0.2);

    float meanPC3dz = hPC3dz->GetBinCenter(hPC3dz->GetMaximumBin());
    float rmsPC3dz = hPC3dz->GetRMS();
    float invPC3dz = hPC3dz->Integral();

    TF1 *fPeakPC3dz = new TF1("fPeakPC3dz", "gaus(0)", -2.5, 3.0);
    fPeakPC3dz->SetParameters(invPC3dz / 50, meanPC3dz, rmsPC3dz);
    fPeakPC3dz->SetLineColor(2);

    hPC3dz->Fit("fPeakPC3dz", "MRQ", "", -2.5, 3.0);

    TPaveText *tPeakPC3dz = new TPaveText(0.435, 0.774, 0.886, 0.848, "brNDC");
    tPeakPC3dz->SetBorderSize(0);
    tPeakPC3dz->SetTextFont(2);
    tPeakPC3dz->SetTextSize(0.04);
    tPeakPC3dz->SetTextColor(1);
    tPeakPC3dz->SetFillColor(10);
    tPeakPC3dz->AddText(Form("Mean=%.3e, Sigma=%.3e", fPeakPC3dz->GetParameter(1), fPeakPC3dz->GetParameter(2)));

    TCanvas *cPc3 = new TCanvas("cPc3", "", 1100, 400);
    cPc3->Divide(2, 1);
    cPc3->cd(1);
    hPC3dphi->Draw();
    tPeakPC3dphi->Draw("SAME");
    species->Draw("SAME");

    cPc3->cd(2);
    hPC3dz->Draw();
    tPeakPC3dz->Draw("SAME");
    cPc3->Update();
    cPc3->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/pAuJet/plots/Pc3.png");
}












