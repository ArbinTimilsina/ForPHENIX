#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/constants.C"

const int NPTBINS = 16;
const float PTBINS[NPTBINS + 1] =
    {
	9.5273, 10.4465, 11.4543, 12.5594, 13.7711, 15.0998, 16.5566, 18.1539, 19.9054,
	21.8258, 23.9315, 26.2404, 28.7720, 31.5479, 34.5915, 37.9289, 41.5882
    };

using namespace std;

void plot510GeV()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TPaveText *speciesCommon = getSpecies(4);
    speciesCommon->AddText("PYTHIA p+p @ 510 GeV");

    TFile *fileParton = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/Run13Simulation/JetHadCorrection/510GeV/parton/JetHadCorrection.root", "r");
    TH1F *hEventsParton = (TH1F*)fileParton->Get("hEvents");
    float nEventsParton = hEventsParton->GetBinContent(1);
    cout << "Parton events: " << nEventsParton << endl;

    TFile *fileParticle = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/Run13Simulation/JetHadCorrection/510GeV/particle/JetHadCorrection.root", "r");
    TH1F *hEventsParticle = (TH1F*)fileParticle->Get("hEvents");
    float nEventsParticle = hEventsParticle->GetBinContent(1);
    cout << "Particle events: " << nEventsParticle << endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH1F *hPtPartonJetsR2 = (TH1F*)fileParton->Get("hPtPartonJetsR2");
    hPtPartonJetsR2->Sumw2();
    hPtPartonJetsR2->Scale(1/nEventsParton, "width");

    TH1F *hPtParticleJetsR2 = (TH1F*)fileParticle->Get("hPtParticleJetsR2");
    hPtParticleJetsR2->Sumw2();
    hPtParticleJetsR2->Scale(1/nEventsParticle, "width");

    TH1F *hRatioR2 = new TH1F("hRatioR2", "", NPTBINS, PTBINS);
    hRatioR2->Sumw2();
    hRatioR2->Divide(hPtPartonJetsR2, hPtParticleJetsR2, 1, 1, "b(1,1) mode");

    TH1F *hPtPartonJetsR3 = (TH1F*)fileParton->Get("hPtPartonJetsR3");
    hPtPartonJetsR3->Sumw2();
    hPtPartonJetsR3->Scale(1/nEventsParton, "width");

    TH1F *hPtParticleJetsR3 = (TH1F*)fileParticle->Get("hPtParticleJetsR3");
    hPtParticleJetsR3->Sumw2();
    hPtParticleJetsR3->Scale(1/nEventsParticle, "width");

    TH1F *hRatioR3 = new TH1F("hRatioR3", "", NPTBINS, PTBINS);
    hRatioR3->Sumw2();
    hRatioR3->Divide(hPtPartonJetsR3, hPtParticleJetsR3, 1, 1, "b(1,1) mode");

    setHisto(hRatioR2, "p_{T} (GeV/c)", "Hadronization Correction Factor");
    hRatioR2->SetMaximum(5.0);
    hRatioR2->SetMinimum(1.0);

    setMarker(hRatioR2, 22, 2, 1.5);
    setMarker(hRatioR3, 23, 4, 1.5);

    TLegend *lRatio = getLegend(3);
    lRatio->AddEntry(hRatioR2, "#color[2]{R = 0.2}", "");
    lRatio->AddEntry(hRatioR3, "#color[4]{R = 0.3}", "");

    TCanvas *cRatio = new TCanvas("cRatio", "", 1100, 700);
    cRatio->cd();
    hRatioR2->Draw("P");
    hRatioR3->Draw("P SAME");

    speciesCommon->Draw("SAME");
    lRatio->Draw("SAME");

    cRatio->Update();
    cRatio->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/plots/JetHadCorrection/510GeV.png");
}































