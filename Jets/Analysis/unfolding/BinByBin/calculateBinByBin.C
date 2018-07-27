#include <TROOT.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>

using namespace std;

void calculateBinByBin(TH2F *hResponseMatrix, TGraphErrors* &gBinByBinEfficiency, TGraphErrors* &gBinByBinPurity, TGraphErrors* &gBinByBinCorrection, bool verbo = false)
{
    cout << endl << endl;
    cout << "*****************************************************************************************************************" << endl;
    cout << "Starting Bin-by-bin calculation. Efficiency, purity, and correction factor will be calculated." << endl << endl;

    const int N = hResponseMatrix->GetNbinsX();
    const int nBins = N;
    cout << "Bins in x-axis and y-axis: " << nBins << endl << endl;

    double sumX[nBins];
    double sumY[nBins];

    double nEfficiencyX[nBins];
    double nEfficiencyErrorX[nBins];
    double nEfficiencyY[nBins];
    double nEfficiencyErrorY[nBins];

    double nPurityX[nBins];
    double nPurityErrorX[nBins];
    double nPurityY[nBins];
    double nPurityErrorY[nBins];

    double nCorrectionX[nBins];
    double nCorrectionErrorX[nBins];
    double nCorrectionY[nBins];
    double nCorrectionErrorY[nBins];

    for(int i = 0; i < nBins; i++)
        {
            sumX[i] = 0.0;
            sumY[i] = 0.0;

            nEfficiencyX[i] = 0.0;
            nEfficiencyErrorX[i] = 0.0;
            nEfficiencyY[i] = 0.0;
            nEfficiencyErrorY[i] = 0.0;

            nPurityX[i] = 0.0;
            nPurityErrorX[i] = 0.0;
            nPurityY[i] = 0.0;
            nPurityErrorY[i] = 0.0;

            nCorrectionX[i] = 0.0;
            nCorrectionErrorX[i] = 0.0;
            nCorrectionY[i] = 0.0;
            nCorrectionErrorY[i] = 0.0;
        }

    for(int ix = 0; ix < nBins; ix++)
        {
            double sum = 0.0;
            for(int iy = 0; iy <= nBins; iy++)
                {
                    //This gives the spectrum of True jet (that were matched), p_{T, true}
                    sum += hResponseMatrix->GetBinContent(ix + 1, iy);
                }
            sumY[ix] = sum;
        }

    for(int iy = 0; iy < nBins; iy++)
        {
            double sum = 0.0;
            for(int ix = 0; ix <= nBins; ix++)
                {
                    //This gives the spectrum of Reco jet (that were matched), p_{T, reco}
                    sum += hResponseMatrix->GetBinContent(ix, iy + 1);
                }
            sumX[iy] = sum;
        }


    for(int ix = 0; ix < nBins; ix++)
        {

            double sameXY = hResponseMatrix->GetBinContent(ix + 1, ix + 1);
            nEfficiencyX[ix] = hResponseMatrix->GetXaxis()->GetBinCenter(ix + 1);
            if(sumY[ix] != 0)
                {
                    nEfficiencyY[ix] = sameXY / sumY[ix];
                    nEfficiencyErrorY[ix] = sqrt(((1 - nEfficiencyY[ix]) * nEfficiencyY[ix]) / sumY[ix]);
                }
            else
                {
                    cout << "sumY is 0 for bin " << ix << endl;
                }

            nPurityX[ix] = hResponseMatrix->GetXaxis()->GetBinCenter(ix + 1);
            if(sumX[ix] != 0)
                {
                    nPurityY[ix] = sameXY / sumX[ix];
                    nPurityErrorY[ix] = sqrt(((1 - nPurityY[ix]) * nPurityY[ix]) / sumX[ix]);
                }
            else
                {
                    cout << "sumX is 0 for bin " << ix << endl;
                }

            nCorrectionX[ix] = hResponseMatrix->GetXaxis()->GetBinCenter(ix + 1);
            if(nEfficiencyY[ix] != 0)
                {
                    nCorrectionY[ix] = nPurityY[ix] / nEfficiencyY[ix];

                    if(sumY[ix] != 0)
                        {
                            nCorrectionErrorY[ix] = (sumX[ix] / pow(sumY[ix], 3)) * (sumX[ix] + sumY[ix] - (2 * sameXY));
                        }
                }
            else
                {
                    cout << "Efficiency is 0 for bin " << ix << endl;
                }
        }

    if(verbo)
        {
            cout << "********************************" << endl;
            cout << "Printing Bin-by-bin factors" << endl;
            cout << "********************************" << endl;
            cout << setw(15) << "Bin Center: " << setw(15) << "E: " << setw(15) << "E E: " << setw(15) << "P: " << setw(15) <<
		"P E: " << setw(15) << "CF: " << setw(15) << "CF E: " << endl;

            for(int ix = 0; ix < nBins; ix++)
                {
                    cout << setw(15) << nEfficiencyX[ix] << setw(15) << nEfficiencyY[ix] << setw(15) << nEfficiencyErrorY[ix] << setw(15) << nPurityY[ix] << setw(15) <<
			nPurityErrorY[ix] << setw(15) << nCorrectionY[ix] << setw(15) << nCorrectionErrorY[ix] << setw(15) << endl;
                }
            cout << "********************************" << endl;
            cout << endl;
        }
    else
        {
            cout << "FYI: Set verbo to true to get bin-by-bin information!!" << endl << endl;
        }

    gBinByBinEfficiency = new TGraphErrors(nBins, nEfficiencyX, nEfficiencyY, nEfficiencyErrorX, nEfficiencyErrorY);
    gBinByBinEfficiency->SetMarkerStyle(33);
    gBinByBinEfficiency->SetMarkerColor(9);

    gBinByBinPurity = new TGraphErrors(nBins, nPurityX, nPurityY, nPurityErrorX, nPurityErrorY);
    gBinByBinPurity->SetMarkerStyle(33);
    gBinByBinPurity->SetMarkerColor(8);

    gBinByBinCorrection = new TGraphErrors(nBins, nCorrectionX, nCorrectionY, nCorrectionErrorX, nCorrectionErrorY);
    gBinByBinCorrection->SetMarkerStyle(33);
    gBinByBinCorrection->SetMarkerColor(2);

    cout << "Bin-by-bin calculations done. Enjoy your graphs!!" << endl;
    cout << "*****************************************************************************************************************" << endl << endl;
}











