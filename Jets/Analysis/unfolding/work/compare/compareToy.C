#include "/direct/phenix+u/arbint/Jets/Analysis/unfolding/BinByBin/doBinByBinUnfolding.C"

#include "/direct/phenix+u/arbint/treasures/unfolding.C"
#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/constants.C"

using namespace std;

void compareToy()
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

    if(false){
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
    }

    ////////////////////////////////////////////////////////////////////////
    //Start here
    ////////////////////////////////////////////////////////////////////////
    TGraphErrors *g1;
    TGraphErrors *g2;
    TGraphErrors *g3;
    TGraphErrors *g4;
    TGraphErrors *g5;

    TGraphErrors *gMyBinByBinUnfolded;
    doBinByBinUnfolding(hMeasured, hMeasuredTrue, hTrue, hResponseMatrix,
                        g1, g2, g3,
                        g4,
                        g5, gMyBinByBinUnfolded);
    gMyBinByBinUnfolded->SetMarkerStyle(33);
    gMyBinByBinUnfolded->SetMarkerColor(46);

    //Correct the Response Matrix for RooUnfold- Truth in y-axis and Reco in x-axis
    flipResponseMatrix(hResponseMatrix);

    TH1F *hRooUnfoldBinByBinUnfolded = getRooUnfoldBinByBin(hMeasured, hMeasuredTrue, hTrue, hResponseMatrix);
    TGraphErrors *gRooUnfoldBinByBinUnfolded = setHistoToGraph(hRooUnfoldBinByBinUnfolded);
    gRooUnfoldBinByBinUnfolded->SetMarkerStyle(33);
    gRooUnfoldBinByBinUnfolded->SetMarkerColor(9);
 
    TH1F *hRooUnfoldBayesUnfolded = getRooUnfoldBayes(hMeasured, hMeasuredTrue, hTrue, hResponseMatrix);
    TGraphErrors *gRooUnfoldBayesUnfolded = setHistoToGraph(hRooUnfoldBayesUnfolded);
    gRooUnfoldBayesUnfolded->SetMarkerStyle(33);
    gRooUnfoldBayesUnfolded->SetMarkerColor(8);

    TH1F *hRooUnfoldSvdUnfolded = getRooUnfoldSvd(hMeasured, hMeasuredTrue, hTrue, hResponseMatrix);
    TGraphErrors *gRooUnfoldSvdUnfolded = setHistoToGraph(hRooUnfoldSvdUnfolded);
    gRooUnfoldSvdUnfolded->SetMarkerStyle(33);
    gRooUnfoldSvdUnfolded->SetMarkerColor(7);

    //hTrue->Sumw2();
    TGraphErrors *gTrue = setHistoToGraph(hTrue);
    gTrue->SetMarkerStyle(33);
    gTrue->SetMarkerColor(2);

    //hMeasured->Sumw2();
    TGraphErrors *gMeasured = setHistoToGraph(hMeasured);
    gMeasured->SetMarkerStyle(33);
    gMeasured->SetMarkerColor(1);

    ////////////////////////////////////////////////////////////////////////
    TMultiGraph *multiGraphCorrected = new TMultiGraph();
    multiGraphCorrected->Add(gTrue);
    multiGraphCorrected->Add(gMeasured);
    multiGraphCorrected->Add(gMyBinByBinUnfolded);
    multiGraphCorrected->Add(gRooUnfoldBinByBinUnfolded);
    multiGraphCorrected->Add(gRooUnfoldBayesUnfolded);
    multiGraphCorrected->Add(gRooUnfoldSvdUnfolded);

    TPaveText *titleCorrected = getHistoTitle();
    titleCorrected->AddText("Comparision of unfolding methods, toy model");

    TLegend *legendCorrected = getLegend(4);
    legendCorrected->AddEntry(gMeasured, "#color[1]{Measured spectra}", " ");
    legendCorrected->AddEntry(gTrue, "#color[2]{True spectra}", " ");
    legendCorrected->AddEntry(gMyBinByBinUnfolded, "#color[46]{My BinByBin Corrected spectra}", " ");
    legendCorrected->AddEntry(gRooUnfoldBinByBinUnfolded, "#color[9]{RooUnfold BinByBin Corrected spectra}", " ");
    legendCorrected->AddEntry(gRooUnfoldBinByBinUnfolded, "#color[8]{RooUnfold Bayes Corrected spectra}", " ");
    legendCorrected->AddEntry(gRooUnfoldBinByBinUnfolded, "#color[7]{RooUnfold SVD Corrected spectra}", " ");

    TCanvas *cCorrected = new TCanvas("cCorrected", "", 1100, 700);
    cCorrected->cd();
    gPad->SetLogy();
    multiGraphCorrected->Draw("AP");
    setMultiGraph(multiGraphCorrected, "p_{T} (GeV/c)", "Count");
    multiGraphCorrected->GetXaxis()->SetLimits(0, 5);
    multiGraphCorrected->SetMinimum(1);
    legendCorrected->Draw("SAME");
    titleCorrected->Draw("SAME");

    cCorrected->Update();
    //cCorrected->SaveAs("/direct/phenix+hhj/arbint/plots/Unfolding/toyModel/CorrectedVariable.png");
}















