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
#include <TCanvas.h>
#include <TCut.h>
#include <TChain.h>
#include <sstream>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#include "/direct/phenix+u/arbint/Jets/Analysis/unfolding/BinByBin/doBinByBinUnfolding.C"

#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"

using namespace std;

void toyModelSimple()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    //***********************************************************************************************************

    TH1F *hTrue = new TH1F("hTrue", "True spectra; p_{T, true}; #frac{dN}{dp_{T}}", 5, 0, 5);
    hTrue->Fill(0.5, 1040);
    hTrue->Fill(1.5, 565);
    hTrue->Fill(2.5, 343);
    hTrue->Fill(3.5, 209);
    hTrue->Fill(4.5, 148);

    /**
    TH1F *hMeasured = new TH1F("hMeasured", "Reconstructed spectra; p_{T, reco}; #frac{dN}{dp_{T}}", 5, 0, 5);
    hMeasured->Fill(0.5, 155);
    hMeasured->Fill(1.5, 133);
    hMeasured->Fill(2.5, 91);
    hMeasured->Fill(3.5, 54);
    hMeasured->Fill(4.5, 29);
    **/

    TH1F *hMeasured = new TH1F("hMeasured", "Reconstructed spectra; p_{T, reco}; #frac{dN}{dp_{T}}", 5, 0, 5);
    hMeasured->Fill(0.5, 255);
    hMeasured->Fill(1.5, 199);
    hMeasured->Fill(2.5, 120);
    hMeasured->Fill(3.5, 65);
    hMeasured->Fill(4.5, 30);

    ////////////////////////////////////////////////////////////////////////
    TH1F *hMeasuredTrue = new TH1F("hMeasuredTrue", "Reconstructed spectra; p_{T, true}; #frac{dN}{dp_{T}}", 5, 0, 5);
    hMeasuredTrue->Fill(0.5, 104);
    hMeasuredTrue->Fill(1.5, 113);
    hMeasuredTrue->Fill(2.5, 103);
    hMeasuredTrue->Fill(3.5, 83);
    hMeasuredTrue->Fill(4.5, 59);


    ////////////////////////////////////////////////////////////////////////
    TH2F *hResponseMatrix = new TH2F("hResponseMatrix", "Response matrix; p_{T, true}; p_{T, reco}",
                                     5, 0, 5, 5, 0, 5);
    hResponseMatrix->Fill(0.5, 0.5, 78.0);
    hResponseMatrix->Fill(1.5, 0.5, 37.0);
    hResponseMatrix->Fill(2.5, 0.5, 22.0);
    hResponseMatrix->Fill(3.5, 0.5, 12.0);
    hResponseMatrix->Fill(4.5, 0.5, 6.0);

    hResponseMatrix->Fill(0.5, 1.5, 21.0);
    hResponseMatrix->Fill(1.5, 1.5, 58.0);
    hResponseMatrix->Fill(2.5, 1.5, 29.0);
    hResponseMatrix->Fill(3.5, 1.5, 16.0);
    hResponseMatrix->Fill(4.5, 1.5, 9.0);

    hResponseMatrix->Fill(0.5, 2.5, 3.0);
    hResponseMatrix->Fill(1.5, 2.5, 15.0);
    hResponseMatrix->Fill(2.5, 2.5, 40.0);
    hResponseMatrix->Fill(3.5, 2.5, 21.0);
    hResponseMatrix->Fill(4.5, 2.5, 12.0);

    hResponseMatrix->Fill(0.5, 3.5, 1.0);
    hResponseMatrix->Fill(1.5, 3.5, 2.0);
    hResponseMatrix->Fill(2.5, 3.5, 10.0);
    hResponseMatrix->Fill(3.5, 3.5, 27.0);
    hResponseMatrix->Fill(4.5, 3.5, 14.0);

    hResponseMatrix->Fill(0.5, 4.5, 1.0);
    hResponseMatrix->Fill(1.5, 4.5, 1.0);
    hResponseMatrix->Fill(2.5, 4.5, 2.0);
    hResponseMatrix->Fill(3.5, 4.5, 7.0);
    hResponseMatrix->Fill(4.5, 4.5, 18.0);

    setHisto(hResponseMatrix, "p_{T, true}", "p_{reco}");

    TPaveText *titleResponse = getHistoTitle();
    titleResponse->AddText("Response Matrix, toy model");

    TCanvas *cResponseMatrix = new TCanvas("cResponseMatrix", "", 700, 700);
    cResponseMatrix->Clear();
    cResponseMatrix->cd();
    hResponseMatrix->Draw("text col");
    titleResponse->Draw("SAME");
    cResponseMatrix->Update();
    cResponseMatrix->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/ResponseMatrixSimple.png");

    ////////////////////////////////////////////////////////////////////////
    //Start here
    TGraphErrors *gEfficiency;
    TGraphErrors *gPurity;
    TGraphErrors *gCorrection;
    TGraphErrors *gMeasuredEfficiency;
    TGraphErrors *gMeasured;
    TGraphErrors *gCorrected;

    doBinByBinUnfolding(hMeasured, 1, 10, hResponseMatrix, hTrue, hMeasuredTrue,
                        gEfficiency, gPurity, gCorrection,
                        gMeasuredEfficiency,
                        gMeasured, gCorrected);

    TLegend *legend = new TLegend(0.81, 0.73, 0.88, 0.88);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetTextFont(2);
    legend->SetTextColor(1);

    legend->AddEntry(gCorrection, "C_{i}", "ep");
    legend->AddEntry(gPurity, "P_{i}", "ep");
    legend->AddEntry(gEfficiency, "E_{i}", "ep");

    TMultiGraph *multiGraph = new TMultiGraph();
    multiGraph->Add(gEfficiency);
    multiGraph->Add(gPurity);
    multiGraph->Add(gCorrection);

    TPaveText *titleBinByBin = getHistoTitle();
    titleBinByBin->AddText("Bin-by-bin quantities, toy model");

    TCanvas *cBinByBin = new TCanvas("cBinByBin", "", 1100, 700);
    cBinByBin->cd();
    multiGraph->Draw("AP");
    setMultiGraph(multiGraph, "p_{T} (GeV/c)", "E_{i} or P_{i} or C_{i}");
    multiGraph->GetXaxis()->SetLimits(0, 5);
    multiGraph->SetMinimum(0);
    multiGraph->SetMaximum(2.8);
    legend->Draw("SAME");
    titleBinByBin->Draw("SAME");

    cBinByBin->Update();
    cBinByBin->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/BinByBinQuantitiesSimple.png");

    ////////////////////////////////////////////////////////////////////////

    setGraph(gMeasuredEfficiency, "p_{T, true} (GeV/c)", "Efficiency");
    gMeasuredEfficiency->SetMinimum(0);
    gMeasuredEfficiency->SetMaximum(0.7);

    TPaveText *titleRecoEfficiency = getHistoTitle();
    titleRecoEfficiency->AddText("Reconstruction Efficiciency, toy model");

    TPaveText *infoEfficiency = getInfo(true);
    infoEfficiency->AddText("Efficiency = #frac{No. of Matched Jets}{No. of True Jets}");

    TCanvas *cRecoEfficiency = new TCanvas("cRecoEfficiency", "", 1100, 700);
    cRecoEfficiency->cd();
    gMeasuredEfficiency->Draw("AP");
    gMeasuredEfficiency->GetXaxis()->SetLimits(0, 5);
    titleRecoEfficiency->Draw("SAME");
    infoEfficiency->Draw("SAME");

    cRecoEfficiency->Update();
    cRecoEfficiency->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/RecoEfficiencySimple.png");


    ////////////////////////////////////////////////////////////////////////
    hTrue->Sumw2();
    TGraphErrors *gTrue = setHistoToGraph(hTrue);
    gTrue->SetMarkerStyle(33);
    gTrue->SetMarkerColor(4);
    gTrue->SetMarkerSize(1.2);

    TGraphErrors *gShiftedCorrected = shiftGraphRight(gCorrected, 0.05);
    TMultiGraph *multiGraphCorrected = new TMultiGraph();
    multiGraphCorrected->Add(gTrue);
    multiGraphCorrected->Add(gMeasured);
    multiGraphCorrected->Add(gShiftedCorrected);

    TPaveText *titleCorrected = getHistoTitle();
    titleCorrected->AddText("Comparision, toy model");

    TLegend *legendCorrected = getLegend(4);
    legendCorrected->AddEntry(gTrue, "True spectra", "ep");
    legendCorrected->AddEntry(gMeasured, "#splitline{Measured spectra}{(uncorrected)}", "ep");
    legendCorrected->AddEntry(gShiftedCorrected, "#splitline{Corrected spectra}{(shifted right)}", "ep");

    TCanvas *cCorrected = new TCanvas("cCorrected", "", 1100, 700);
    cCorrected->cd();
    gPad->SetLogy();
    multiGraphCorrected->Draw("AP");
    setMultiGraph(multiGraphCorrected, "p_{T} (GeV/c)", "#frac{dN}{dp_{T}}");
    multiGraphCorrected->GetXaxis()->SetLimits(0, 5);
    legendCorrected->Draw("SAME");
    titleCorrected->Draw("SAME");

    cCorrected->Update();
    cCorrected->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/CorrectedSimple.png");
}














