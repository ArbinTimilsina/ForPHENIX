#include "calculateBinByBin.C"
#include "/direct/phenix+u/arbint/treasures/histoMakeup.C"
#include "/direct/phenix+u/arbint/treasures/calculateEfficiency.C"

using namespace std;

void doBinByBinUnfolding(TH1F *hMeasured, TH1F *hTrue, TH2F *hResponseMatrix,
                         TGraphErrors* &gBinByBinEfficiency, TGraphErrors* &gBinByBinPurity, TGraphErrors* &gBinByBinCorrection,
                         TGraphErrors* &gReconstructionEfficiency,
                         TGraphErrors* &gMeasured, TGraphErrors* &gCorrected,
                         bool verbo = false)
{
    cout << endl;
    cout << "*****************************************************************************************************************" << endl;
    cout << "Starting Bin-by-bin Unfolding." << endl;
    cout << "*****************************************************************************************************************" << endl;

    ///////////////////////////////////////////////
    //Get Bin-by-bin correction factors
    ///////////////////////////////////////////////
    calculateBinByBin(hResponseMatrix, gBinByBinEfficiency, gBinByBinPurity, gBinByBinCorrection, verbo);

    ///////////////////////////////////////////////
    //Get Reconstrucion Efficiency
    ///////////////////////////////////////////////
    const int binsTrue = hTrue->GetNbinsX();
    TH1F *hPass = (TH1F*)hTrue->Clone("hPass");
    hPass->Clear();
    for(int ix = 0; ix < binsTrue; ix++)
        {
            float sum = 0.0;
            float sumE = 0.0;
            for(int iy = 0; iy <= binsTrue; iy++)
                {
                    //This gives the spectrum of True jet (that were matched), p_{T, true}
                    sum += hResponseMatrix->GetBinContent(ix + 1, iy);
		    sumE += pow(hResponseMatrix->GetBinError(ix + 1, iy),2);
                }
	    hPass->SetBinContent(ix + 1, sum);
	    hPass->SetBinError(ix + 1, sqrt(sumE));
        }
    calculateEfficiency(hTrue, hPass, gReconstructionEfficiency, verbo);

    ///////////////////////////////////////////////
    //Do the correction
    ///////////////////////////////////////////////
    double *nBinByBinCorrectionY      = gBinByBinCorrection->GetY();
    double *nBinByBinCorrectionErrorY = gBinByBinCorrection->GetEY();

    double *nReconstructionEfficiencyY      = gReconstructionEfficiency->GetY();
    double *nReconstructionEfficiencyErrorY = gReconstructionEfficiency->GetEY();

    gMeasured = setHistoToGraph(hMeasured, false);
    double *nMeasuredY      = gMeasured->GetY();
    double *nMeasuredErrorY = gMeasured->GetEY();

    double *nMeasuredX      = gMeasured->GetX();
    const int N = gMeasured->GetN();
    const int nBins = N;

    double nCorrectedX[nBins];
    double nCorrectedErrorX[nBins];
    double nCorrectedY[nBins];
    double nCorrectedErrorY[nBins];

    for(int i = 0; i < nBins; i++)
        {
            nCorrectedX[i] = 0.0;
            nCorrectedErrorX[i] = 0.0;
            nCorrectedY[i] = 0.0;
            nCorrectedErrorY[i] = 0.0;
        }

    for(int ix = 0; ix < nBins; ix++)
        {
            nCorrectedX[ix] = nMeasuredX[ix];

            if(nReconstructionEfficiencyY[ix] != 0)
                {
                    nCorrectedY[ix] = (1 / nReconstructionEfficiencyY[ix]) * nBinByBinCorrectionY[ix] * nMeasuredY[ix];
                    
                    if((nBinByBinCorrectionY[ix] != 0) && (nMeasuredY[ix] != 0))
                        {
                            nCorrectedErrorY[ix] = fabs(nCorrectedY[ix]) *
				sqrt(pow((nReconstructionEfficiencyErrorY[ix] / nReconstructionEfficiencyY[ix]), 2) +
				     pow((nBinByBinCorrectionErrorY[ix] / nBinByBinCorrectionY[ix]), 2) +
				     pow((nMeasuredErrorY[ix] / nMeasuredY[ix]), 2));
                        }
                    else
                        {
                            cout << "Either Measured or Bin-By-Bin Correction is 0 for bin " << ix << endl;
                        }
                }
            else
                {
                    cout << "Reconstruction Efficiency is 0 for bin " << ix << endl;
                }
        }

    gCorrected = new TGraphErrors(nBins, nCorrectedX, nCorrectedY, nCorrectedErrorX, nCorrectedErrorY);
    gCorrected->SetMarkerStyle(33);
    gCorrected->SetMarkerColor(2);

    gMeasured->SetMarkerStyle(33);
    gMeasured->SetMarkerColor(1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    if(verbo)
        {
            cout << endl;
            cout << "*****************************************************************************************************************" << endl;
            cout << "Bin-by-bin Unfolding done. Enjoy your graphs!!" << endl;
            cout << "*****************************************************************************************************************" << endl;

            cout << "********************************" << endl;
            cout << "Printing Measured information" << endl;
            cout << "********************************" << endl;
            cout << setw(15) << "Bin Center: " << setw(15) << "Measured: " << setw(15) << "Error: " << setw(15)
                 << "Error %: " << setw(15) << endl;
            for(int ix = 0; ix < nBins; ix++)
                {
                    if(nMeasuredY[ix] != 0)
                        {
                            cout << setw(15) << nMeasuredX[ix] << setw(15) << nMeasuredY[ix] << setw(15) << nMeasuredErrorY[ix]
                                 << setw(15) << nMeasuredErrorY[ix] * 100 / nMeasuredY[ix] << setw(15) << endl;
                        }
                    else
                        {
                            cout << setw(15) << nMeasuredX[ix] << setw(15) << nMeasuredY[ix] << setw(15) << nMeasuredErrorY[ix]
                                 << setw(15) << "NA" << setw(15) << endl;
                        }
                }
            cout << "********************************" << endl << endl;


            cout << "********************************" << endl;
            cout << "Printing Corrected information" << endl;
            cout << "********************************" << endl;
            cout << setw(15) << "Bin Center: " << setw(15) << "Corrected: " << setw(15) << "Error: " << setw(15)
                 << "Error %: " << setw(15) << endl;
            for(int ix = 0; ix < nBins; ix++)
                {
                    if(nCorrectedY[ix] != 0)
                        {
                            cout << setw(15) << nCorrectedX[ix] << setw(15) << nCorrectedY[ix] << setw(15) << nCorrectedErrorY[ix]
                                 << setw(15) << nCorrectedErrorY[ix] * 100 / nCorrectedY[ix] << setw(15) << endl;
                        }
                    else
                        {
                            cout << setw(15) << nCorrectedX[ix] << setw(15) << nCorrectedY[ix] << setw(15) << nCorrectedErrorY[ix]
                                 << setw(15) << "NA" << setw(15) << endl;
                        }
                }
            cout << "********************************" << endl;
            cout << "*****************************************************************************************************************" << endl;
            cout << endl;
        }
    else
        {
            cout << "FYI: Set verbo to true to get bin-by-bin information!!" << endl << endl;
        }

}























