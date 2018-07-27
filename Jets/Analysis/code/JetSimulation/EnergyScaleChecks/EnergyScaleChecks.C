//General PHENIX tools
#include "PHCompositeNode.h"
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include "phool.h"
#include <getClass.h>
#include "RunHeader.h"

//Fun4All tools
#include "Fun4AllServer.h"
#include "Fun4AllHistoManager.h"
#include "Fun4AllReturnCodes.h"

//Root tools
#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include <TNtuple.h>

#include <VtxOut.h>

//My source file
#include "EnergyScaleChecks.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
EnergyScaleChecks::EnergyScaleChecks()
    : SubsysReco("EnergyScaleChecks")
{
    nTotalEvents = 0;

    return;
}

int EnergyScaleChecks::Init(PHCompositeNode *topNode)
{
    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(true);
    if(isCuAu)
        {
            jetAnalyzer->SetCuAu(true);
        }
    else
        {
            jetAnalyzer->SetCuAu(false);
        }
    jetAnalyzer->MyInit();

    //Output file name
    outfile = new TFile("EnergyScaleChecks.root", "RECREATE");

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    for (unsigned int s = 0; s < 8; s++)
        {
            hPi0Mass[s] = new TH1F(Form("hPi0Mass_%u", s), "Pi0 mass", 720, 0.05, 0.35);
        }

    //*******************************************************************************************************
    return EVENT_OK;
}

int EnergyScaleChecks::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int EnergyScaleChecks::ResetEvent(PHCompositeNode *topNode)
{
    return EVENT_OK;
}


int EnergyScaleChecks::process_event(PHCompositeNode *topNode)
{
    nTotalEvents++;

    VtxOut *vertexOut = getClass<VtxOut>(topNode, "VtxOut");
    if (!vertexOut)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    float zvertex = vertexOut->get_ZVertex();

    //Cluster
    std::vector<clusters> cluster_list;
    cluster_list.clear();
    jetAnalyzer->SetData(false);
    jetAnalyzer->GetClusters(topNode, cluster_list, zvertex);

    GetPi0(cluster_list, hPi0Mass);

    return 0;
}

void EnergyScaleChecks::SetCuAu(bool what = true)
{
    isCuAu = what;
}

void EnergyScaleChecks::GetPi0(std::vector<clusters> cluster_list,
                               TH1F *hMass[8])
{
    std::vector<photons> all_photons;
    for (unsigned int c = 0; c < cluster_list.size(); c++)
        {
            bool passGoodCluster = cluster_list[c].passIsValid && cluster_list[c].passNotBad;
            bool passEnergy = cluster_list[c].energy > 0.5;

            int armsect = cluster_list[c].armsect;
            bool passChi2 = !EmcMap::IsPbGl(armsect) && cluster_list[c].chi2 < 3.0;
            bool passProb = EmcMap::IsPbGl(armsect) && cluster_list[c].prob > 0.02;

            bool isPhoton = passGoodCluster && passEnergy && (passChi2 || passProb);

            if(isPhoton)
                {
                    photons tempPhoton;
                    tempPhoton.armsect = cluster_list[c].armsect;
                    tempPhoton.pT = cluster_list[c].pT;
                    tempPhoton.eta = cluster_list[c].eta;
                    tempPhoton.phi = cluster_list[c].phi;

                    all_photons.push_back(tempPhoton);
                }
        }

    for (unsigned int p = 0; p < all_photons.size(); p++)
        {
            for (unsigned int pp = p + 1; pp < all_photons.size(); pp++)
                {
                    //Require same armsect
                    if(all_photons[p].armsect != all_photons[pp].armsect)
                        {
                            continue;
                        }

                    float asymmetry = fabs((all_photons[p].pT - all_photons[pp].pT) / (all_photons[p].pT + all_photons[pp].pT));
                    if(asymmetry > 0.8)
                        {
                            continue;
                        }

                    TLorentzVector photon1;
                    photon1.SetPtEtaPhiM(all_photons[p].pT, all_photons[p].eta, all_photons[p].phi, 0.0);
                    TLorentzVector photon2;
                    photon2.SetPtEtaPhiM(all_photons[pp].pT, all_photons[pp].eta, all_photons[pp].phi, 0.0);

                    TLorentzVector pion = photon1 + photon2;

                    float pionPt = pion.Pt();
                    if(pionPt < 3.0 || pionPt > 10.0)
                        {
                            continue;
                        }

                    float pionMass = pion.M();
                    hMass[all_photons[p].armsect]->Fill(pionMass);
                }
        }
}

int EnergyScaleChecks::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int EnergyScaleChecks::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nTotalEvents);

    cout << endl;
    cout << "Total events processed: " << nTotalEvents << endl;
    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}

