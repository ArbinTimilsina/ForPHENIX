#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/constants.C"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "RooUnfoldResponse.h"
#include "RooUnfold.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBayes.h"

//Needed for plotting d-vector
#include "TSVDUnfold.h"
#endif

using namespace std;

void doUnfoldingR3()
{
#ifdef __CINT__
    gSystem->Load("libRooUnfold");
#endif

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TPaveText *algorithm = new TPaveText(0.728, 0.830, 0.897, 0.881, "brNDC");
    algorithm->SetTextSize(0.035);
    algorithm->SetBorderSize(1);
    algorithm->SetTextFont(2);
    algorithm->SetFillColor(10);
    algorithm->SetTextColor(46);
    algorithm->SetTextAlign(32);
    algorithm->AddText("Anti-kt, R = 0.3");

    TPaveText *speciesUnfolded = getSpecies(3);
    speciesUnfolded->SetTextSize(0.05);

    TFile *fileFakeSubtracted = new TFile("/phenix/hhj/arbint/RootFiles/PostQM15/RootFiles/FakeSubtractedR3.root", "READ");
    TH1F *hFakeSubtracted = (TH1F*)fileFakeSubtracted->Get("hFakeSubtractedR3");

    TFile *fileSim = new TFile("/phenix/hhj/arbint/RootFiles/PostQM15/SimR3/JetAnalyzerSimR3.root", "r");

    TH1F *hSimTrue = (TH1F*)fileSim->Get("hPtTrueJets_0");
    TH2F *hResponseMatrix = (TH2F*)fileSim->Get("hResponseMatrixDefault_0");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Do unfolding: SVD
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Error are selected with RooUnfold::ErrorTreatment enum:
    //No error treatment (RooUnfold::kNoError)- (uses sqrt(N))
    //Bin-by-bin errors (RooUnfold::kErrors)- (no correlations)
    //Full covariance matrix (RooUnfold::kCovariance) propagated through the unfolding
    //Covariance matrix from the variation of the results in toy MC tests (RooUnfold::kCovToy).
    //This last method should be more accurate, especially for RooUnfoldBayes.
    RooUnfold::ErrorTreatment errorSVD = RooUnfold::kCovariance;

    //"measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
    //but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
    //in "truth" for unmeasured events (inefficiency).

    //TH1F *hSimTrueNull = (TH1F*)fileSim->Get("hPtTrueJets_0");
    //hSimTrueNull->Reset();

    TH1F *hSimRecoNull = (TH1F*)fileSim->Get("hPtRecoJets_0");
    hSimRecoNull->Reset();

    int kRegSVD = 3;
    RooUnfoldResponse response(hSimRecoNull, hSimTrue, hResponseMatrix);

    RooUnfoldSvd unfoldSvd(&response, hFakeSubtracted, kRegSVD);
    TH1F *hSVDUnfolded = (TH1F*) unfoldSvd.Hreco(errorSVD);

    //For d-vector
    TSVDUnfold *myTSVD = (TSVDUnfold*) unfoldSvd.Impl();
    TH1D *hDVector = myTSVD->GetD();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Plot stuff
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH1F *hSVDUnfoldedClone = (TH1F*)hSVDUnfolded->Clone("hSVDUnfoldedClone");
    setMarkerWithColorError(hFakeSubtracted, 33, 2, 2.0);
    setMarkerWithColorError(hSVDUnfoldedClone, 33, 1, 2.0);

    hSVDUnfoldedClone->GetXaxis()->SetTitle("p_{T} or p_{T, Reco} (GeV/c)");
    hSVDUnfoldedClone->GetXaxis()->SetTitleSize(0.04);
    hSVDUnfoldedClone->GetXaxis()->SetTitleOffset(1.15);
    hSVDUnfoldedClone->GetXaxis()->SetTitleFont(2);
    hSVDUnfoldedClone->GetXaxis()->SetTickLength(0.04);
    hSVDUnfoldedClone->GetXaxis()->SetLabelFont(2);

    hSVDUnfoldedClone->GetYaxis()->SetTitle("Per-event Counts/Bin");
    hSVDUnfoldedClone->GetYaxis()->SetTitleSize(0.05);
    hSVDUnfoldedClone->GetYaxis()->SetTitleOffset(0.8);
    hSVDUnfoldedClone->GetYaxis()->SetTitleFont(2);
    hSVDUnfoldedClone->GetYaxis()->SetTickLength(0.04);
    hSVDUnfoldedClone->GetYaxis()->SetLabelFont(2);
    hSVDUnfoldedClone->GetYaxis()->SetLabelSize(0.04);
    hSVDUnfoldedClone->SetTitle("");

    //hSVDUnfoldedClone->SetMaximum(histoMax);
    hSVDUnfoldedClone->SetMinimum(2E-9);

    TLegend *legendMethod = getLegend(3);
    legendMethod->AddEntry(hFakeSubtracted, "#color[2]{\"Fake\" Jets Subtracted}", "");
    legendMethod->AddEntry(hSVDUnfoldedClone, Form("#color[1]{SVD unfolded (kReg = %d)}", kRegSVD), "");

    TCanvas *cUnfoldingMethod = new TCanvas("cUnfoldingMethos", "", 1100, 700);
    cUnfoldingMethod->cd();

    gPad->SetLogy();
    hSVDUnfoldedClone->Draw("P");
    hFakeSubtracted->Draw("P SAME");
    legendMethod->Draw("SAME");
    //speciesUnfolded->Draw("SAME");
    algorithm->Draw("SAME");
    cUnfoldingMethod->Update();
    cUnfoldingMethod->SaveAs("/phenix/hhj/arbint/plots/PostQM15/R3/unfolding/Unfolded.png");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Plot d-Vector
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    TLine *lineVectorX = new TLine(0, 1, 19, 1);
    lineVectorX->SetLineStyle(2);
    lineVectorX->SetLineColor(1);
    lineVectorX->SetLineWidth(2);

    TLine *lineVectorY = new TLine(kRegSVD, 0, kRegSVD, 0.8 * hDVector->GetMaximum());
    lineVectorY->SetLineStyle(1);
    lineVectorY->SetLineColor(2);
    lineVectorY->SetLineWidth(2);

    setHisto(hDVector, "i", "|d_{i}|");

    hDVector->GetXaxis()->SetTitleSize(0.045);
    hDVector->GetYaxis()->SetTitleSize(0.045);

    hDVector->SetLineWidth(2);
    hDVector->SetMinimum(5E-2);

    TPaveText *textVector = new TPaveText(0.328, 0.731, 0.498, 0.783, "brNDC");
    textVector->SetTextSize(0.045);
    textVector->SetBorderSize(0);
    textVector->SetTextFont(2);
    textVector->SetFillColor(10);
    textVector->SetTextColor(2);
    textVector->SetTextAlign(12);
    textVector->AddText(Form("kReg = %d", kRegSVD));

    TCanvas *cCanvasVector = new TCanvas("cCanvasVector", "", 1100, 700);
    cCanvasVector->cd();
    gPad->SetLogy();
    hDVector->Draw();
    lineVectorX->Draw("SAME");
    lineVectorY->Draw("SAME");
    textVector->Draw("SAME");
    TPaveText *title = getHistoTitle();
    title->AddText("R=0.3");
    title->Draw("SAME");
    cCanvasVector->SaveAs("/phenix/hhj/arbint/plots/PostQM15/R3/unfolding/dVector.png");

    TFile *fileSave = new TFile("/phenix/hhj/arbint/RootFiles/PostQM15/RootFiles/UnfoldedR3.root", "RECREATE");
    hSVDUnfolded->SetName("hSVDUnfoldedR3");
    hSVDUnfolded->SetTitle("SVD unfolded, R=0.3");
    hSVDUnfolded->Write();
    fileSave->Write();
    fileSave->Close();
}





















