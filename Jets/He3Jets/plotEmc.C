#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void plotEmc()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/He3Jet.root", "READ");;
    TPaveText *species = getSpecies(3);
    species->AddText("Run 14 He3+Au @ 200 GeV, ERT");

    //***************************************************************************************
    //***************************************************************************************

    TH1F *hEMCdphi = (TH1F*)file->Get("hEMCdphi");
    hEMCdphi->Rebin(8);
    hEMCdphi->GetXaxis()->SetRangeUser(-0.04, 0.04);
    hEMCdphi->SetLineColor(1);
    hEMCdphi->SetLineWidth(2);
    setHisto(hEMCdphi, "emcdphi", "Counts/Bin");
    setHistoMax(hEMCdphi, false, 0.2);

    float meanEMCdphi = hEMCdphi->GetBinCenter(hEMCdphi->GetMaximumBin());
    float rmsEMCdphi = hEMCdphi->GetRMS();
    float invEMCdphi = hEMCdphi->Integral();

    TF1 *fPeakEMCdphi = new TF1("fPeakEMCdphi", "gaus(0)", -0.015, 0.015);
    fPeakEMCdphi->SetParameters(invEMCdphi / 50, meanEMCdphi, rmsEMCdphi);
    fPeakEMCdphi->SetLineColor(2);

    hEMCdphi->Fit("fPeakEMCdphi", "MRQ", "", -0.015, 0.015);

    TPaveText *tPeakEMCdphi = new TPaveText(0.435, 0.774, 0.886, 0.848, "brNDC");
    tPeakEMCdphi->SetBorderSize(0);
    tPeakEMCdphi->SetTextFont(2);
    tPeakEMCdphi->SetTextSize(0.04);
    tPeakEMCdphi->SetTextColor(1);
    tPeakEMCdphi->SetFillColor(10);
    tPeakEMCdphi->AddText(Form("Mean=%.3e, Sigma=%.3e", fPeakEMCdphi->GetParameter(1), fPeakEMCdphi->GetParameter(2)));

    TH1F *hEMCdz = (TH1F*)file->Get("hEMCdz");
    hEMCdz->GetXaxis()->SetRangeUser(-15.0, 15.0);
    hEMCdz->SetLineColor(1);
    hEMCdz->SetLineWidth(2);
    setHisto(hEMCdz, "emcdz", "Counts/Bin");
    setHistoMax(hEMCdz, false, 0.2);

    float meanEMCdz = hEMCdz->GetBinCenter(hEMCdz->GetMaximumBin());
    float rmsEMCdz = hEMCdz->GetRMS();
    float invEMCdz = hEMCdz->Integral();

    TF1 *fPeakEMCdz = new TF1("fPeakEMCdz", "gaus(0)", -5.0, 5.0);
    fPeakEMCdz->SetParameters(invEMCdz / 50, meanEMCdz, rmsEMCdz);
    fPeakEMCdz->SetLineColor(2);

    hEMCdz->Fit("fPeakEMCdz", "MRQ", "", -5.0, 5.0);

    TPaveText *tPeakEMCdz = new TPaveText(0.435, 0.774, 0.886, 0.848, "brNDC");
    tPeakEMCdz->SetBorderSize(0);
    tPeakEMCdz->SetTextFont(2);
    tPeakEMCdz->SetTextSize(0.04);
    tPeakEMCdz->SetTextColor(1);
    tPeakEMCdz->SetFillColor(10);
    tPeakEMCdz->AddText(Form("Mean=%.3e, Sigma=%.3e", fPeakEMCdz->GetParameter(1), fPeakEMCdz->GetParameter(2)));

    TCanvas *cEmc = new TCanvas("cEmc", "", 1100, 400);
    cEmc->Divide(2, 1);
    cEmc->cd(1);
    hEMCdphi->Draw();
    tPeakEMCdphi->Draw("SAME");
    species->Draw("SAME");

    cEmc->cd(2);
    hEMCdz->Draw();
    tPeakEMCdz->Draw("SAME");
    cEmc->Update();
    cEmc->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/plots/Emc.png");
}












