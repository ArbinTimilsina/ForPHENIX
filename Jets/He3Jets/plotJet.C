#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void plotJet()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/He3Jet.root", "READ");;
    TPaveText *species = getSpecies(3);
    species->AddText("He3+Au @ 200 GeV, anti-kt, R = 0.3 jet");

    TH1F *hCf = (TH1F*)file->Get("hCf");
    hCf->SetLineColor(1);
    hCf->SetLineWidth(2);
    setHisto(hCf, "Charged fraction", "Counts");
    setHistoMax(hCf, false, 0.2);

    TPaveText *tCf= new TPaveText(0.435, 0.774, 0.886, 0.848, "brNDC");
    tCf->SetBorderSize(0);
    tCf->SetTextFont(2);
    tCf->SetTextSize(0.04);
    tCf->SetTextColor(1);
    tCf->SetFillColor(10);
    tCf->AddText("Jet-level cuts: pT > 6.0 + n.c. > = 3");

    TCanvas *cCf = new TCanvas("cCf", "", 1100, 700);
    cCf->Clear();
    cCf->cd();
    gPad->SetLogy();
    hCf->Draw();
    tCf->Draw("SAME");
    species->Draw("SAME");
    cCf->Update();
    cCf->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/plots/Cf.png");


    TNtuple *tJets = (TNtuple*)file->Get("tJets");

    TH1F *hJetNoCut = new TH1F("hJetNoCut", "", 60, 6.0, 66.0);
    hJetNoCut->SetLineColor(1);
    hJetNoCut->SetLineWidth(2);
    setHisto(hJetNoCut, "p_{T, Reco} (GeV/c)", "Counts/bin");
    hJetNoCut->SetMinimum(1.0);

    TH1F *hJetCut1 = new TH1F("hJetCut1", "", 60, 6.0, 66.0);
    hJetCut1->SetLineColor(4);
    hJetCut1->SetLineWidth(2);

    TH1F *hJetCut2 = new TH1F("hJetCut2", "", 60, 6.0, 66.0);
    hJetCut2->SetLineColor(2);
    hJetCut2->SetLineWidth(2);

    tJets->Draw("pT>>hJetNoCut");
    tJets->Draw("pT>>hJetCut1" ,"cf > 0.1 && cf < 0.9");
    tJets->Draw("pT>>hJetCut2", "cf > 0.2 && cf < 0.8");

    TLegend *lJets = getLegend(3);
    lJets->AddEntry(hJetNoCut, "#color[1]{nc>=3}", "");
    lJets->AddEntry(hJetCut1, "#color[4]{nc>=3 && cf>0.1 && cf<0.9}", "");
    lJets->AddEntry(hJetCut2, "#color[2]{nc>=3 && cf>0.2 && cf<0.8}", "");

    TCanvas *cJets = new TCanvas("cJets", "", 1100, 700);
    cJets->cd();
    gPad->SetLogy();

    hJetNoCut->Draw();
    hJetCut1->Draw("SAME");
    hJetCut2->Draw("SAME");
    lJets->Draw("SAME");
    species->Draw("SAME");

    cJets->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/plots/Jets.png");
}

