#include "/direct/phenix+u/arbint/Jets/Analysis/unfolding/BinByBin/doBinByBinUnfolding.C"

#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/constants.C"

#include "/direct/phenix+u/arbint/treasures/unfolding.C"

using namespace std;

void comparePP()
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
    TH1F *hReco = (TH1F*)fileSim->Get("hPtPisaRecoJet_0");
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

    doBinByBinUnfolding(hPtMB, hReco, hTrue, hResponseMatrix,
                        g1, g2, g3,
                        g4,
                        gMeasured, gMyBinByBinUnfolded,
			true);

    TCanvas *c = new TCanvas("c", " ", 1100, 700);
    c->cd();
    g3->Draw("APE");
    //g4->SetMaximum(0.65);
    //g4->SetMinimum(0.0);

    gMeasured->SetMarkerStyle(33);
    gMeasured->SetMarkerColor(1);

    gMyBinByBinUnfolded->SetMarkerStyle(33);
    gMyBinByBinUnfolded->SetMarkerColor(2);

    //RooUnfold
    flipResponseMatrix(hResponseMatrix);
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
    TGraphErrors *gRooUnfoldBinByBinUnfoldedShifted = shiftGraphRight(gRooUnfoldBinByBinUnfolded,0.3);
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
    legendCompare->AddEntry(gRooUnfoldBinByBinUnfoldedShifted, "#color[4]{RooUnfold Bin-By-Bin Unfolded (shifted right)}", "");
    legendCompare->AddEntry(gRooUnfoldBayesUnfoldedShifted, "#color[8]{RooUnfold Bayes' Unfolded (1 iteration)}", "");

    TMultiGraph *multiGraphCompare = new TMultiGraph();
    multiGraphCompare->Add(gMeasured);
    multiGraphCompare->Add(gMyBinByBinUnfoldedShifted);
    multiGraphCompare->Add(gRooUnfoldBinByBinUnfoldedShifted);
    multiGraphCompare->Add(gRooUnfoldBayesUnfoldedShifted);

    TPaveText *titleComparision = getHistoTitle();
    titleComparision->AddText("Comparision of Unfolding Methods");

    cCompare->cd();
    gPad->SetLogy();
    multiGraphCompare->Draw("AP");
    legendCompare->Draw("SAME");
    titleComparision->Draw("SAME");
    setMultiGraph(multiGraphCompare, "p_{T} (GeV/c)", "Yield");
    multiGraphCompare->GetXaxis()->SetLimits(MINPT, MAXPT);
    multiGraphCompare->SetMinimum(0.5);
    cCompare->Update();
}





