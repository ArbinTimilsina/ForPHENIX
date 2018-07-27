#include "/direct/phenix+u/arbint/Jets/Analysis/unfolding/BinByBin/doBinByBinUnfolding.C"

#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/constants.C"

#include "/direct/phenix+u/arbint/treasures/unfolding.C"

using namespace std;

TH2F* getNormalizedMatrix(TH1F *hMeasured, TH2F *hResponseMatrix)
{
    const int N = hResponseMatrix->GetNbinsX();
    const int nBins = N;

    hMeasured->Sumw2();
    float nSpectraDataX[nBins];
    float nSpectraDataErrorX[nBins];
    float nSpectraDataY[nBins];
    float nSpectraDataErrorY[nBins];
    for(int ix = 0; ix < nBins; ix++)
        {
            nSpectraDataX[ix] = hMeasured->GetXaxis()->GetBinCenter(ix + 1);
            nSpectraDataErrorX[ix] = 0.0;
            nSpectraDataY[ix] = hMeasured->GetBinContent(ix + 1);
            nSpectraDataErrorY[ix] = hMeasured->GetBinError(ix + 1);
        }


    int firstBin = hMeasured->GetXaxis()->FindBin(20.0);
    int lastBin  = hMeasured->GetXaxis()->FindBin(26.0);//26.5);

    float nScaleData = hMeasured->Integral(firstBin - 1, lastBin - 1);
    float nScaleMC   = 0.0;
    for(int iy = 1; iy <= nBins; iy++)
        {
            for(int ix = firstBin; ix <= lastBin; ix++)
                {
                    nScaleMC += hResponseMatrix->GetBinContent(ix, iy);
                }
        }


    float nSpectraMCX[nBins];
    float nSpectraMCErrorX[nBins];
    float nSpectraMCY[nBins];
    float nSpectraMCErrorY[nBins];
    for(int iy = 0; iy < nBins; iy++)
        {
            float sumX = 0.0;
            float sumEX = 0.0;
            for(int ix = 0; ix < nBins; ix++)
                {
                    sumX += hResponseMatrix->GetBinContent(ix + 1, iy + 1);
                    sumEX += pow(hResponseMatrix->GetBinError(ix + 1, iy + 1), 2);
                }

            nSpectraMCX[iy] = hResponseMatrix->GetXaxis()->GetBinCenter(iy + 1);
            nSpectraMCErrorX[iy] = 0.0;
            nSpectraMCY[iy] = sumX * (nScaleData / nScaleMC);
            nSpectraMCErrorY[iy] = sqrt(sumEX * (nScaleData / nScaleMC));
        }


    float nRatioX[nBins];
    float nRatioErrorX[nBins];
    float nRatioY[nBins];
    float nRatioErrorY[nBins];
    for(int ix = 0; ix < nBins; ix++)
        {
            nRatioX[ix] = nSpectraDataX[ix];
            nRatioErrorX[ix] = 0.0;
            if(nSpectraMCY[ix] != 0 && nSpectraDataY[ix])
                {
                    nRatioY[ix] = nSpectraDataY[ix] / nSpectraMCY[ix];
                    nRatioErrorY[ix] = fabs(nRatioY[ix]) *
			sqrt(pow((nSpectraDataErrorY[ix] / nSpectraDataY[ix]), 2) +
			     pow((nSpectraMCErrorY[ix] / nSpectraMCY[ix]), 2));
                }
            else
                {
                    nRatioY[ix] = 1.0;
                    nRatioErrorY[ix] = 0.0;
                }
        }

    TH2F *hResponseMatrixNormalized = (TH2F*)hResponseMatrix->Clone("hResponseMatrixNormalized");
    for(int iy = 0; iy < nBins; iy++)
        {
            for(int ix = 0; ix < nBins; ix++)
                {
                    double normalized = hResponseMatrix->GetBinContent(ix + 1, iy + 1) * nRatioY[iy + 1];
                    hResponseMatrixNormalized->SetBinContent(ix + 1, iy + 1, normalized);
                }
        }

    float nSpectraMCCorrectedX[nBins];
    float nSpectraMCCorrectedErrorX[nBins];
    float nSpectraMCCorrectedY[nBins];
    float nSpectraMCCorrectedErrorY[nBins];
    for(int iy = 0; iy < nBins; iy++)
        {
            float sumX = 0.0;
            float sumEX = 0.0;
            for(int ix = 0; ix < nBins; ix++)
                {
                    sumX += hResponseMatrixNormalized->GetBinContent(ix+1, iy+1);
                    sumEX += pow(hResponseMatrixNormalized->GetBinError(ix+1, iy+1), 2);
                }

            nSpectraMCCorrectedX[iy] = hResponseMatrixNormalized->GetXaxis()->GetBinCenter(iy+1);
            nSpectraMCCorrectedErrorX[iy] = 0.0;
	    nSpectraMCCorrectedY[iy] = sumX * (nScaleData / nScaleMC);
            nSpectraMCCorrectedErrorY[iy] = sqrt(sumEX * (nScaleData / nScaleMC));
        }


    //plot stuff
    if(false){
    gDataYield = new TGraphErrors(nBins, nSpectraDataX, nSpectraDataY, nSpectraDataErrorX, nSpectraDataErrorY);
    gMCYield = new TGraphErrors(nBins, nSpectraMCX, nSpectraMCY, nSpectraMCErrorX, nSpectraMCErrorY);
    gRatio = new TGraphErrors(nBins, nRatioX, nRatioY, nRatioErrorX, nRatioErrorY);
    gMCYieldCorrected = new TGraphErrors(nBins, nSpectraMCCorrectedX, nSpectraMCCorrectedY, nSpectraMCCorrectedErrorX, nSpectraMCCorrectedErrorY);

    gDataYield->SetMarkerStyle(33);
    gDataYield->SetMarkerColor(1);

    gMCYield->SetMarkerStyle(33);
    gMCYield->SetMarkerColor(2);

    gMCYieldCorrected->SetMarkerStyle(33);
    gMCYieldCorrected->SetMarkerColor(3);

    TMultiGraph *multiGraph1 = new TMultiGraph();
    multiGraph1->Add(gDataYield);
    multiGraph1->Add(gMCYield);
    multiGraph1->Add(gMCYieldCorrected);

    TLegend *legend1 = getLegend(3);
    legend1->AddEntry(gDataYield, "#color[1]{Data}", "");
    legend1->AddEntry(gMCYield, "#color[2]{MC}", "");
    legend1->AddEntry(gMCYieldCorrected, "#color[3]{MC Corrected}", "");

    TCanvas *cCanvas1 = new TCanvas("cCanvas1", "", 1100, 700);
    cCanvas1->cd();
    gPad->SetLogy();
    multiGraph1->Draw("APZ");
    setMultiGraph(multiGraph1, "p_{T} (GeV/c)", "Yield");
    //multiGraph1->GetXaxis()->SetLimits(MINPT, MAXPT);
    //multiGraph1->SetMaximum(1.35);
    multiGraph1->SetMinimum(0.1);

    legend1->Draw("SAME");
    cCanvas1->Update();
    }

    return hResponseMatrixNormalized;
}


void compareMethods()
{
#ifdef __CINT__
    gSystem->Load("libRooUnfold");
#endif

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);


    //***************************************************************************************************
    //Get the p+p plots
    //***************************************************************************************************
    TFile *file = new TFile("/direct/phenix+hhj/arbint/RootFiles/JetAnalyzerPP_ERT.root", "READ");
    TH1F *hPtMB = (TH1F*)file->Get("hPtCent0");

    //***************************************************************************************************
    //Get hisots needed for unfolding
    //***************************************************************************************************
    TFile *fileSim = new TFile("/direct/phenix+hhj/arbint/JetSimulation/JetAnalyzerSim/PP/JetAnalyzerSim.root", "r");
    TH2F *hResponseMatrix = (TH2F*)fileSim->Get("hPisaRecoResponseMatrix_0");
    TH1F *hTrue = (TH1F*)fileSim->Get("hPtTrueJet_0");
    TH1F *hReco = (TH1F*)fileSim->Get("hPtPisaRecoJe_0");
    TH1F *hPass = (TH1F*)fileSim->Get("hForPisaRecoMatchingEfficiency_0");

    //TH2F *hResponseMatrix = getNormalizedMatrix(hPtMB, hResponseMatrix1);

    //***************************************************************************************************
    //Do the unfolding
    //***************************************************************************************************
    //My Method
    TGraphErrors *gMeasured;
    TGraphErrors *gMyBinByBinUnfolded;

    TGraphErrors *g1;
    TGraphErrors *g2;
    TGraphErrors *g3;
    TGraphErrors *g4;

    doBinByBinUnfolding(hPtMB, hResponseMatrix, hTrue, hPass,
                        g1, g2, g3,
                        g4,
                        gMeasured, gMyBinByBinUnfolded);

    gMeasured->SetMarkerStyle(33);
    gMeasured->SetMarkerColor(1);

    gMyBinByBinUnfolded->SetMarkerStyle(33);
    gMyBinByBinUnfolded->SetMarkerColor(2);

    flipResponseMatrix(hResponseMatrix);

    //RooUnfold Bin-By-Bin
    TH1F *hRooUnfoldBinByBinUnfolded = getRooUnfoldBinByBin(hPtMB, hReco, hTrue, hResponseMatrix);
    TGraphErrors *gRooUnfoldBinByBinUnfolded = setHistoToGraph(hRooUnfoldBinByBinUnfolded);
    gRooUnfoldBinByBinUnfolded->SetMarkerStyle(33);
    gRooUnfoldBinByBinUnfolded->SetMarkerColor(4);

    TH1F *hRooUnfoldBayesUnfolded = getRooUnfoldBayes(hPtMB, hReco, hTrue, hResponseMatrix);
    TGraphErrors *gRooUnfoldBayesUnfolded = setHistoToGraph(hRooUnfoldBayesUnfolded);
    gRooUnfoldBayesUnfolded->SetMarkerStyle(33);
    gRooUnfoldBayesUnfolded->SetMarkerColor(8);

    //***************************************************************************************************
    //Compare
    //***************************************************************************************************
    TCanvas *cCompare = new TCanvas("cCompare", " ", 1100, 700);
    cCompare->Clear();

    TGraphErrors *gMyBinByBinUnfoldedShifted = shiftGraphRight(gMyBinByBinUnfolded,0.0);
    TGraphErrors *gRooUnfoldBinByBinUnfoldedShifted = shiftGraphRight(gRooUnfoldBinByBinUnfolded,0.0);
    TGraphErrors *gRooUnfoldBayesUnfoldedShifted = shiftGraphRight(gRooUnfoldBayesUnfolded,0.0);

    TLegend *legendCompare = new TLegend(0.62, 0.66, 0.90, 0.81);
    legendCompare->SetTextSize(0.03);
    legendCompare->SetBorderSize(0);
    legendCompare->SetFillColor(10);
    legendCompare->SetTextFont(2);
    legendCompare->SetTextColor(1);
    legendCompare->SetTextAlign(32);
    legendCompare->AddEntry(gMeasured, "#color[1]{Not corrected p+p}", "");
    legendCompare->AddEntry(gMyBinByBinUnfoldedShifted, "#color[2]{My Bin-By-Bin Unfolded}", "");
    legendCompare->AddEntry(gRooUnfoldBinByBinUnfoldedShifted, "#color[4]{RooUnfold Bin-By-Bin Unfolded}", "");
    legendCompare->AddEntry(gRooUnfoldBayesUnfoldedShifted, "#color[8]{RooUnfold Bayes' Unfolded (1 iteration)}", "");

    TMultiGraph *multiGraphCompare = new TMultiGraph();
    multiGraphCompare->Add(gMeasured);
    multiGraphCompare->Add(gMyBinByBinUnfoldedShifted);
    multiGraphCompare->Add(gRooUnfoldBinByBinUnfoldedShifted);
    multiGraphCompare->Add(gRooUnfoldBayesUnfoldedShifted);

    TPaveText *titleComparision = getHistoTitle();
    titleComparision->AddText("Comparision of Unfolding Methods- inital look");

    cCompare->cd();
    gPad->SetLogy();
    multiGraphCompare->Draw("AP");
    legendCompare->Draw("SAME");
    titleComparision->Draw("SAME");
    setMultiGraph(multiGraphCompare, "p_{T} (GeV/c)", "Yield");
    multiGraphCompare->GetXaxis()->SetLimits(10.0, 50.0);
    multiGraphCompare->SetMinimum(0.8);
    cCompare->Update();
}





