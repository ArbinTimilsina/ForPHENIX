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

void toyModelVariable()
{
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    //***********************************************************************************************************
    const int pTBins = 11;
    const float pT[pTBins + 1] = {0.0, 0.2, 0.4, 0.6, 0.9, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

    TH1F *hTrue = new TH1F("hTrue", "True spectra; p_{T, true}; Count", pTBins, pT);
    hTrue->Fill(0.1, 6465);
    hTrue->Fill(0.3, 4894);
    hTrue->Fill(0.5, 3533);
    hTrue->Fill(0.7, 2463);
    hTrue->Fill(1.0, 1659);
    hTrue->Fill(1.3, 1090);
    hTrue->Fill(1.6, 649);
    hTrue->Fill(2.1, 428);
    hTrue->Fill(2.6, 257);
    hTrue->Fill(3.5, 147);
    hTrue->Fill(4.5, 82);


    TH1F *hMeasured = new TH1F("hMeasured", "Reconstructed spectra; p_{T, reco}; #frac{dN}{dp_{T}}", pTBins, pT);
    hMeasured->Fill(0.1, 601);
    hMeasured->Fill(0.3, 522);
    hMeasured->Fill(0.5, 383);
    hMeasured->Fill(0.7, 265);
    hMeasured->Fill(1.0, 177);
    hMeasured->Fill(1.3, 114);
    hMeasured->Fill(1.6, 72);
    hMeasured->Fill(2.1, 45);
    hMeasured->Fill(2.6, 29);
    hMeasured->Fill(3.5, 18);
    hMeasured->Fill(4.5, 13);

    ////////////////////////////////////////////////////////////////////////
    TH1F *hMeasuredTrue = new TH1F("hMeasuredTrue", "Reconstructed spectra; p_{T, true}; #frac{dN}{dp_{T}}", pTBins, pT);
    hMeasuredTrue->Fill(0.1, 386);
    hMeasuredTrue->Fill(0.3, 425);
    hMeasuredTrue->Fill(0.5, 386);
    hMeasuredTrue->Fill(0.7, 316);
    hMeasuredTrue->Fill(1.0, 242);
    hMeasuredTrue->Fill(1.3, 175);
    hMeasuredTrue->Fill(1.6, 122);
    hMeasuredTrue->Fill(2.1, 81);
    hMeasuredTrue->Fill(2.6, 52);
    hMeasuredTrue->Fill(3.5, 33);
    hMeasuredTrue->Fill(4.5, 21);

    ////////////////////////////////////////////////////////////////////////
    TH2F *hResponseMatrix = new TH2F("hResponseMatrix", "Response matrix; p_{T, true}; p_{T, reco}",
                                     pTBins, pT, pTBins, pT);
    hResponseMatrix->Fill(0.1, 0.1, 286);
    hResponseMatrix->Fill(0.3, 0.1, 142);
    hResponseMatrix->Fill(0.5, 0.1, 81);
    hResponseMatrix->Fill(0.7, 0.1, 44);
    hResponseMatrix->Fill(1.0, 0.1, 23);
    hResponseMatrix->Fill(1.3, 0.1, 12);
    hResponseMatrix->Fill(1.6, 0.1, 6);
    hResponseMatrix->Fill(2.1, 0.1, 3);
    hResponseMatrix->Fill(2.6, 0.1, 2);
    hResponseMatrix->Fill(3.5, 0.1, 1);
    hResponseMatrix->Fill(4.5, 0.1, 1);

    hResponseMatrix->Fill(0.1, 0.3, 76);
    hResponseMatrix->Fill(0.3, 0.3, 211);
    hResponseMatrix->Fill(0.5, 0.3, 107);
    hResponseMatrix->Fill(0.7, 0.3, 60);
    hResponseMatrix->Fill(1.0, 0.3, 33);
    hResponseMatrix->Fill(1.3, 0.3, 17);
    hResponseMatrix->Fill(1.6, 0.3, 9);
    hResponseMatrix->Fill(2.1, 0.3, 5);
    hResponseMatrix->Fill(2.6, 0.3, 2);
    hResponseMatrix->Fill(3.5, 0.3, 1);
    hResponseMatrix->Fill(4.5, 0.3, 1);

    hResponseMatrix->Fill(0.1, 0.5, 13);
    hResponseMatrix->Fill(0.3, 0.5, 54);
    hResponseMatrix->Fill(0.5, 0.5, 148);
    hResponseMatrix->Fill(0.7, 0.5, 78);
    hResponseMatrix->Fill(1.0, 0.5, 43);
    hResponseMatrix->Fill(1.3, 0.5, 23);
    hResponseMatrix->Fill(1.6, 0.5, 12);
    hResponseMatrix->Fill(2.1, 0.5, 6);
    hResponseMatrix->Fill(2.6, 0.5, 3);
    hResponseMatrix->Fill(3.5, 0.5, 2);
    hResponseMatrix->Fill(4.5, 0.5, 1);

    hResponseMatrix->Fill(0.1, 0.7, 4);
    hResponseMatrix->Fill(0.3, 0.7, 9);
    hResponseMatrix->Fill(0.5, 0.7, 37);
    hResponseMatrix->Fill(0.7, 0.7, 100);
    hResponseMatrix->Fill(1.0, 0.7, 54);
    hResponseMatrix->Fill(1.3, 0.7, 30);
    hResponseMatrix->Fill(1.6, 0.7, 16);
    hResponseMatrix->Fill(2.1, 0.7, 8);
    hResponseMatrix->Fill(2.6, 0.7, 4);
    hResponseMatrix->Fill(3.5, 0.7, 2);
    hResponseMatrix->Fill(4.5, 0.7, 1);

    hResponseMatrix->Fill(0.1, 1.0, 1);
    hResponseMatrix->Fill(0.3, 1.0, 3);
    hResponseMatrix->Fill(0.5, 1.0, 6);
    hResponseMatrix->Fill(0.7, 1.0, 25);
    hResponseMatrix->Fill(1.0, 1.0, 66);
    hResponseMatrix->Fill(1.3, 1.0, 37);
    hResponseMatrix->Fill(1.6, 1.0, 20);
    hResponseMatrix->Fill(2.1, 1.0, 10);
    hResponseMatrix->Fill(2.6, 1.0, 5);
    hResponseMatrix->Fill(3.5, 1.0, 3);
    hResponseMatrix->Fill(4.5, 1.0, 1);

    hResponseMatrix->Fill(0.1, 1.3, 1);
    hResponseMatrix->Fill(0.3, 1.3, 1);
    hResponseMatrix->Fill(0.5, 1.3, 2);
    hResponseMatrix->Fill(0.7, 1.3, 4);
    hResponseMatrix->Fill(1.0, 1.3, 16);
    hResponseMatrix->Fill(1.3, 1.3, 41);
    hResponseMatrix->Fill(1.6, 1.3, 24);
    hResponseMatrix->Fill(2.1, 1.3, 13);
    hResponseMatrix->Fill(2.6, 1.3, 7);
    hResponseMatrix->Fill(3.5, 1.3, 3);
    hResponseMatrix->Fill(4.5, 1.3, 2);

    hResponseMatrix->Fill(0.1, 1.6, 1);
    hResponseMatrix->Fill(0.3, 1.6, 1);
    hResponseMatrix->Fill(0.5, 1.6, 1);
    hResponseMatrix->Fill(0.7, 1.6, 1);
    hResponseMatrix->Fill(1.0, 1.6, 3);
    hResponseMatrix->Fill(1.3, 1.6, 10);
    hResponseMatrix->Fill(1.6, 1.6, 26);
    hResponseMatrix->Fill(2.1, 1.6, 15);
    hResponseMatrix->Fill(2.6, 1.6, 8);
    hResponseMatrix->Fill(3.5, 1.6, 4);
    hResponseMatrix->Fill(4.5, 1.6, 2);

    hResponseMatrix->Fill(0.1, 2.1, 1);
    hResponseMatrix->Fill(0.3, 2.1, 1);
    hResponseMatrix->Fill(0.5, 2.1, 1);
    hResponseMatrix->Fill(0.7, 2.1, 1);
    hResponseMatrix->Fill(1.0, 2.1, 1);
    hResponseMatrix->Fill(1.3, 2.1, 2);
    hResponseMatrix->Fill(1.6, 2.1, 6);
    hResponseMatrix->Fill(2.1, 2.1, 15);
    hResponseMatrix->Fill(2.6, 2.1, 9);
    hResponseMatrix->Fill(3.5, 2.1, 5);
    hResponseMatrix->Fill(4.5, 2.1, 3);

    hResponseMatrix->Fill(0.1, 2.6, 1);
    hResponseMatrix->Fill(0.3, 2.6, 1);
    hResponseMatrix->Fill(0.5, 2.6, 1);
    hResponseMatrix->Fill(0.7, 2.6, 1);
    hResponseMatrix->Fill(1.0, 2.6, 1);
    hResponseMatrix->Fill(1.3, 2.6, 1);
    hResponseMatrix->Fill(1.6, 2.6, 1);
    hResponseMatrix->Fill(2.1, 2.6, 4);
    hResponseMatrix->Fill(2.6, 2.6, 9);
    hResponseMatrix->Fill(3.5, 2.6, 6);
    hResponseMatrix->Fill(4.5, 2.6, 3);

    hResponseMatrix->Fill(0.1, 3.5, 1);
    hResponseMatrix->Fill(0.3, 3.5, 1);
    hResponseMatrix->Fill(0.5, 3.5, 1);
    hResponseMatrix->Fill(0.7, 3.5, 1);
    hResponseMatrix->Fill(1.0, 3.5, 1);
    hResponseMatrix->Fill(1.3, 3.5, 1);
    hResponseMatrix->Fill(1.6, 3.5, 1);
    hResponseMatrix->Fill(2.1, 3.5, 1);
    hResponseMatrix->Fill(2.6, 3.5, 2);
    hResponseMatrix->Fill(3.5, 3.5, 5);
    hResponseMatrix->Fill(4.5, 3.5, 3);

    hResponseMatrix->Fill(0.1, 4.5, 1);
    hResponseMatrix->Fill(0.3, 4.5, 1);
    hResponseMatrix->Fill(0.5, 4.5, 1);
    hResponseMatrix->Fill(0.7, 4.5, 1);
    hResponseMatrix->Fill(1.0, 4.5, 1);
    hResponseMatrix->Fill(1.3, 4.5, 1);
    hResponseMatrix->Fill(1.6, 4.5, 1);
    hResponseMatrix->Fill(2.1, 4.5, 1);
    hResponseMatrix->Fill(2.6, 4.5, 1);
    hResponseMatrix->Fill(3.5, 4.5, 1);
    hResponseMatrix->Fill(4.5, 4.5, 3);

    setHisto(hResponseMatrix, "p_{T, true}", "p_{T, reco}");
    hResponseMatrix->GetXaxis()->SetTickLength(0.01);
    hResponseMatrix->GetYaxis()->SetTickLength(0.01);

    TPaveText *titleResponse = getHistoTitle();
    titleResponse->AddText("Response Matrix, toy model");

    TCanvas *cResponseMatrix = new TCanvas("cResponseMatrix", "", 750, 750);
    cResponseMatrix->Clear();
    cResponseMatrix->cd();
    gPad->SetLogz();
    hResponseMatrix->Draw("text col");
    titleResponse->Draw("SAME");
    cResponseMatrix->Update();
    cResponseMatrix->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/ResponseMatrixVariable.png");

    ////////////////////////////////////////////////////////////////////////
    //Start here
    TGraphErrors *gEfficiency;
    TGraphErrors *gPurity;
    TGraphErrors *gCorrection;
    TGraphErrors *gMeasuredEfficiency;
    TGraphErrors *gMeasured;
    TGraphErrors *gCorrected;

    doBinByBinUnfolding(hMeasured, hResponseMatrix, hTrue, hMeasuredTrue,
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
    cBinByBin->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/BinByBinQuantitiesVariable.png");

    ////////////////////////////////////////////////////////////////////////

    setGraph(gMeasuredEfficiency, "p_{T, true} (GeV/c)", "Efficiency");
    gMeasuredEfficiency->SetMinimum(0);
    gMeasuredEfficiency->SetMaximum(0.5);

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
    cRecoEfficiency->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/RecoEfficiencyVariable.png");


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
    setMultiGraph(multiGraphCorrected, "p_{T} (GeV/c)", "Count");
    multiGraphCorrected->GetXaxis()->SetLimits(0, 5);
    legendCorrected->Draw("SAME");
    titleCorrected->Draw("SAME");

    cCorrected->Update();
    cCorrected->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/CorrectedVariable.png");
}















