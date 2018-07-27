#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void plotPt()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file = new TFile("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/He3Jet.root", "READ");;
    TPaveText *species = getSpecies(3);
    species->AddText("He3+Au @ 200 GeV, anti-kt, R = 0.3 jet");

    TNtuple *tJets = (TNtuple*)file->Get("tJets");

    TH1F *hJet = new TH1F("hJet", "", 60, 6.0, 66.0);
    hJet->SetLineColor(1);
    hJet->SetLineWidth(2);
    setHisto(hJet, "p_{T} (GeV/c)", "Counts/bin");
    hJet->SetMaximum(1E6);
    hJet->SetMinimum(0.5);

    tJets->Draw("pT>>hJet", "cf > 0.2 && cf < 0.8");

    TLegend *lJets = getLegend(3);
    lJets->AddEntry(hJet, "#color[2]{nc>=3 && cf>0.2 && cf<0.8}", "");

    TCanvas *cJets = new TCanvas("cJets", "", 1100, 700);
    cJets->cd();
    gPad->SetLogy();
    hJet->Draw();
    species->Draw("SAME");
    lJets->Draw("SAME");

    cJets->SaveAs("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/arbint/He3Jet/plots/CountVsPt.png");

    cout<<endl<<endl;
    cout << setw(15) << "Bin center" << setw(15) << "count" << endl;
    for(int ix = 1; ix <= 60; ix++)
        {
            cout << setw(15) << hJet->GetBinCenter(ix) << setw(15) << hJet->GetBinContent(ix) << endl;
        }
}

