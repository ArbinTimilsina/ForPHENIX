//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <PHObject.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TNtuple.h>

//Pythia tools
#include <PHHijingHeader.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//My source file
#include "JetAnalyzerSim.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
JetAnalyzerSim::JetAnalyzerSim(const float ir, const float iminPt,
                               const float inc, const float iminCf, const float imaxCf,
                               const float iminDeltaR)
    : SubsysReco("JetAnalyzerSim")
{
    R               = ir;
    minPt           = iminPt;

    nc              = inc;
    minCf           = iminCf;
    maxCf           = imaxCf;

    minDeltaR       = iminDeltaR;

    nTotalEvents = 0;
    nTotalEventsCentrality1 = 0;
    nTotalEventsCentrality2 = 0;
    nTotalEventsCentrality3 = 0;
    nTotalEventsCentrality4 = 0;

    treesFill         = false;

    return;
}

int JetAnalyzerSim::Init(PHCompositeNode *topNode)
{
    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");

    //Output file name
    outfile = new TFile("JetAnalyzerSim.root", "RECREATE");

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    hEvents = new TH1F("hEvents", "Number of events", 50, 0, 50);

    for (unsigned int c = 0; c < 5; c++)
        {
            hPtTrueJets[c] = new TH1F(Form("hPtTrueJets_%u", c), "p_{T} of True Reconstructed Jets of PYTHIA",
                                      NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtTrueJetsEast[c] = new TH1F(Form("hPtTrueJetsEast_%u", c), "p_{T} of True Reconstructed Jets of PYTHIA (East arm)",
                                          NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtTrueJetsWest[c] = new TH1F(Form("hPtTrueJetsWest_%u", c), "p_{T} of True Reconstructed Jets of PYTHIA (West arm)",
                                          NPTBINS_TRUE_R3, PTBINS_TRUE_R3);


            hPtRecoJets[c] = new TH1F(Form("hPtRecoJets_%u", c), "p_{T} of Reco Jets",
                                      NPTBINS_RECO_R3, PTBINS_RECO_R3);


            hDistance[c] = new TH1F(Form("hDistance_%u", c), "#DeltaR between True Jets and Reco/Embedded Jets",
                                    400, 0.0, 0.4);
            hDistanceVsPtTrue[c] = new TH2F(Form("hDistanceVsPtTrue_%u", c), "#DeltaR vs p_{T, true} for Reco/Embedded Jets",
                                            40, 10, 50, 40, 0.0, 0.4);
            hPtTrueMatchedJets[c] = new TH1F(Form("hPtTrueMatchedJets_%u", c), "p_{T} of True Matched Jets",
                                             NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedJets[c] = new TH1F(Form("hPtRecoMatchedJets_%u", c), "p_{T} of Reco Matched Jets",
                                             NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrix[c] = new TH2F(Form("hResponseMatrix_%u", c), "p_{T, true} vs p_{T, reco} for Reco/Embedded Jets",
                                          NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            for (unsigned int r = 0; r < 22; r++)
                {
                    hForJESR[c][r] = new TH1F(Form("hForJESR_%u_%u", c, r), "For JES and JER",
                                              120, 0.0, 4.0);
                }
        }

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Variations
    //*******************************************************************************************************
    for (unsigned int c = 0; c < 5; c++)
        {
            hPtRecoDefault[c] = new TH1F(Form("hPtRecoDefault_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedDefault[c] = new TH1F(Form("hPtTrueMatchedDefault_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedDefault[c] = new TH1F(Form("hPtRecoMatchedDefault_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixDefault[c] = new TH2F(Form("hResponseMatrixDefault_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoEast[c] = new TH1F(Form("hPtRecoEast_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedEast[c] = new TH1F(Form("hPtTrueMatchedEast_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedEast[c] = new TH1F(Form("hPtRecoMatchedEast_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixEast[c] = new TH2F(Form("hResponseMatrixEast_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoWest[c] = new TH1F(Form("hPtRecoWest_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedWest[c] = new TH1F(Form("hPtTrueMatchedWest_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedWest[c] = new TH1F(Form("hPtRecoMatchedWest_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixWest[c] = new TH2F(Form("hResponseMatrixWest_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoFidTight[c] = new TH1F(Form("hPtRecoFidTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedFidTight[c] = new TH1F(Form("hPtTrueMatchedFidTight_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedFidTight[c] = new TH1F(Form("hPtRecoMatchedFidTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixFidTight[c] = new TH2F(Form("hResponseMatrixFidTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoNc[c] = new TH1F(Form("hPtRecoNc_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedNc[c] = new TH1F(Form("hPtTrueMatchedNc_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedNc[c] = new TH1F(Form("hPtRecoMatchedNc_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixNc[c] = new TH2F(Form("hResponseMatrixNc_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoCf[c] = new TH1F(Form("hPtRecoCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedCf[c] = new TH1F(Form("hPtTrueMatchedCf_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedCf[c] = new TH1F(Form("hPtRecoMatchedCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixCf[c] = new TH2F(Form("hResponseMatrixCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoNcCf[c] = new TH1F(Form("hPtRecoNcCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedNcCf[c] = new TH1F(Form("hPtTrueMatchedNcCf_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedNcCf[c] = new TH1F(Form("hPtRecoMatchedNcCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixNcCf[c] = new TH2F(Form("hResponseMatrixNcCf_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoTrackPlus[c] = new TH1F(Form("hPtRecoTrackPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedTrackPlus[c] = new TH1F(Form("hPtTrueMatchedTrackPlus_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedTrackPlus[c] = new TH1F(Form("hPtRecoMatchedTrackPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixTrackPlus[c] = new TH2F(Form("hResponseMatrixTrackPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoTrackMinus[c] = new TH1F(Form("hPtRecoTrackMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedTrackMinus[c] = new TH1F(Form("hPtTrueMatchedTrackMinus_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedTrackMinus[c] = new TH1F(Form("hPtRecoMatchedTrackMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixTrackMinus[c] = new TH2F(Form("hResponseMatrixTrackMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoClusterPlus[c] = new TH1F(Form("hPtRecoClusterPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedClusterPlus[c] = new TH1F(Form("hPtTrueMatchedClusterPlus_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedClusterPlus[c] = new TH1F(Form("hPtRecoMatchedClusterPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixClusterPlus[c] = new TH2F(Form("hResponseMatrixClusterPlus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoClusterMinus[c] = new TH1F(Form("hPtRecoClusterMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedClusterMinus[c] = new TH1F(Form("hPtTrueMatchedClusterMinus_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedClusterMinus[c] = new TH1F(Form("hPtRecoMatchedClusterMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixClusterMinus[c] = new TH2F(Form("hResponseMatrixClusterMinus_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);

            hPtRecoTrClTight[c] = new TH1F(Form("hPtRecoTrClTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hPtTrueMatchedTrClTight[c] = new TH1F(Form("hPtTrueMatchedTrClTight_%u", c), "", NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
            hPtRecoMatchedTrClTight[c] = new TH1F(Form("hPtRecoMatchedTrClTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hResponseMatrixTrClTight[c] = new TH2F(Form("hResponseMatrixTrClTight_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3, NPTBINS_TRUE_R3, PTBINS_TRUE_R3);
        }

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    if(treesFill)
        {
            tTrueJets = new TTree("tTrueJets", "True Jets of Pythia");
            tTrueJets->Branch("tTrueJets", &true_jets);

            //*******************************************************************************************************

            tTotalParticles = new TTree("tTotalParticles", "Total Tracks and Clusters for Jet Reconstruction");
            tTotalParticles->Branch("tTotalParticles", &total_particles);

            //*******************************************************************************************************

            tRecoJets = new TTree("tRecoJets", "Reconstrctred Jets");
            tRecoJets->Branch("tRecoJets", &reco_jets);

            tMatchedJets = new TTree("tMatchedJets", "Matched Jets");
            tMatchedJets->Branch("tMatchedJets", &matched_jets);


            //*******************************************************************************************************
        }

    return EVENT_OK;
}

int JetAnalyzerSim::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int JetAnalyzerSim::ResetEvent(PHCompositeNode *topNode)
{
    true_jets.clear();

    total_particles.clear();

    reco_jets.clear();
    matched_jets.clear();

    return EVENT_OK;
}


int JetAnalyzerSim::process_event(PHCompositeNode *topNode)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Event info
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    PHHijingHeader *eventContainer = findNode::getClass<PHHijingHeader>(topNode, "EventContainerNode");
    if (!eventContainer)
        {
            cout << "No EventContainer! No sense continuing" << endl;
            exit(1);
        }
    float centrality = eventContainer->GetBimpact();
    centralityBin = jetAnalyzer->getCentralityBin(centrality);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //True jet info
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    PHPythiaContainer *jetContainer = findNode::getClass<PHPythiaContainer>(topNode, "JetContainerNode");
    if (!jetContainer)
        {
            cout << "No JetContainer! No sense continuing" << endl;
            exit(1);
        }

    for (int iJet = 0; iJet < jetContainer->size(); iJet++)
        {
            TMCParticle *jet = jetContainer->getParticle(iJet);

            jets tempTrueJet;
            tempTrueJet.arm           = jet->GetKF();
            tempTrueJet.centralityBin = centralityBin;
            tempTrueJet.pT            = jet->GetPx();
            tempTrueJet.eta           = jet->GetPy();
            tempTrueJet.phi           = jet->GetPz();
            tempTrueJet.nc            = jet->GetVx();
            tempTrueJet.cf            = jet->GetVy();
            tempTrueJet.nf            = jet->GetVz();
            tempTrueJet.disc          = jet->GetTime();

            true_jets.push_back(tempTrueJet);
        }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Particle info
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    PHPythiaContainer *particleContainer = findNode::getClass<PHPythiaContainer>(topNode, "ParticleContainerNode");
    if (!particleContainer)
        {
            cout << "No ParticleContainer! No sense continuing" << endl;
            exit(1);
        }

    for (int iParticle = 0; iParticle < particleContainer->size(); iParticle++)
        {
            TMCParticle *particle = particleContainer->getParticle(iParticle);

            particles tempParticle;
            tempParticle.arm          = particle->GetKF();
            tempParticle.charge       = particle->GetKS();
            tempParticle.energy       = particle->GetEnergy();
            tempParticle.mom          = particle->GetVx();
            tempParticle.pT           = particle->GetVy();
            tempParticle.eT           = particle->GetVz();
            tempParticle.px           = particle->GetPx();
            tempParticle.py           = particle->GetPy();
            tempParticle.pz           = particle->GetPz();
            tempParticle.eta          = particle->GetMass();
            tempParticle.phi          = particle->GetTime();
            tempParticle.phiDC        = particle->GetLifetime();

            tempParticle.ertTrigger   = false;
            tempParticle.id           = 0;
            tempParticle.jetPt        = -9999.9;

            total_particles.push_back(tempParticle);
        }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    nTotalEvents++;
    if(centralityBin == 1)
        {
            nTotalEventsCentrality1++;
        }
    if(centralityBin == 2)
        {
            nTotalEventsCentrality2++;
        }
    if(centralityBin == 3)
        {
            nTotalEventsCentrality3++;
        }
    if(centralityBin == 4)
        {
            nTotalEventsCentrality4++;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Get True Jet
    ///////////////////////////////////////////////////////////////////////////////////////////////
    getPythiaTrueJet(true_jets,
                     hPtTrueJets, hPtTrueJetsEast, hPtTrueJetsWest);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    getRecoJet(total_particles, reco_jets, hPtRecoJets);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    getMatchedJet(minDeltaR,
                  true_jets, reco_jets, matched_jets,
                  hDistance, hDistanceVsPtTrue, hPtTrueMatchedJets, hPtRecoMatchedJets, hResponseMatrix, hForJESR);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Variations
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////
    //Default
    //////////////////////////////////////////////////
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, total_particles,
                 hPtRecoDefault, hPtTrueMatchedDefault, hPtRecoMatchedDefault, hResponseMatrixDefault);

    //////////////////////////////////////////////////
    //East arm
    //////////////////////////////////////////////////
    getVariation(nc, minCf, maxCf, 1,
                 true_jets, total_particles,
                 hPtRecoEast, hPtTrueMatchedEast, hPtRecoMatchedEast, hResponseMatrixEast);

    //////////////////////////////////////////////////
    //West arm
    //////////////////////////////////////////////////
    getVariation(nc, minCf, maxCf, 2,
                 true_jets, total_particles,
                 hPtRecoWest, hPtTrueMatchedWest, hPtRecoMatchedWest, hResponseMatrixWest);

    //////////////////////////////////////////////////
    //Tight fiducial cut
    //////////////////////////////////////////////////
    getVariation(nc, minCf, maxCf, 3,
                 true_jets, total_particles,
                 hPtRecoFidTight, hPtTrueMatchedFidTight, hPtRecoMatchedFidTight, hResponseMatrixFidTight);

    //////////////////////////////////////////////////
    //nc>=5
    //////////////////////////////////////////////////
    getVariation(5.0, minCf, maxCf, 0,
                 true_jets, total_particles,
                 hPtRecoNc, hPtTrueMatchedNc, hPtRecoMatchedNc, hResponseMatrixNc);

    //////////////////////////////////////////////////
    //cf<0.6
    //////////////////////////////////////////////////
    getVariation(nc, minCf, 0.6, 0,
                 true_jets, total_particles,
                 hPtRecoCf, hPtTrueMatchedCf, hPtRecoMatchedCf, hResponseMatrixCf);

    //////////////////////////////////////////////////
    //nc>=5 && cf<0.6
    //////////////////////////////////////////////////
    getVariation(5.0, minCf, 0.6, 0,
                 true_jets, total_particles,
                 hPtRecoNcCf, hPtTrueMatchedNcCf, hPtRecoMatchedNcCf, hResponseMatrixNcCf);

    //////////////////////////////////////////////////
    //Tracks pT varied by +- 2% < 10 GeV/c...
    //////////////////////////////////////////////////

    //////////////////////////////////////////////////
    //Plus
    std::vector<particles> track_plus;
    track_plus.clear();
    for(int i = 0; i < total_particles.size(); i++)
        {
            particles tempParticle;
            tempParticle.arm          = total_particles[i].arm;
            tempParticle.charge       = total_particles[i].charge;
            tempParticle.energy       = total_particles[i].energy;
            tempParticle.eT           = total_particles[i].eT;
            tempParticle.eta          = total_particles[i].eta;
            tempParticle.phi          = total_particles[i].phi;
            tempParticle.phiDC        = total_particles[i].phiDC;
            tempParticle.ertTrigger   = total_particles[i].ertTrigger;
            tempParticle.id           = total_particles[i].id;
            tempParticle.jetPt        = total_particles[i].jetPt;

            if(tempParticle.charge != 0)
                {
                    float factor = 0.02;
                    if(total_particles[i].pT > 10)
                        {
                            factor = 0.02 + 0.01 * ((total_particles[i].pT - 10.0) / 10.0);
                        }
                    tempParticle.pT           = total_particles[i].pT + (factor * total_particles[i].pT);
                    tempParticle.mom          = tempParticle.pT / (total_particles[i].pT / total_particles[i].mom);
                    tempParticle.px           = tempParticle.pT * cos(total_particles[i].phi);
                    tempParticle.py           = tempParticle.pT * sin(total_particles[i].phi);
                    tempParticle.pz           = total_particles[i].pz;
                }
            else
                {
                    tempParticle.pT           = total_particles[i].pT;
                    tempParticle.mom          = total_particles[i].mom;
                    tempParticle.px           = total_particles[i].px;
                    tempParticle.py           = total_particles[i].py;
                    tempParticle.pz           = total_particles[i].pz;
                }
            track_plus.push_back(tempParticle);
        }
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, track_plus,
                 hPtRecoTrackPlus, hPtTrueMatchedTrackPlus, hPtRecoMatchedTrackPlus, hResponseMatrixTrackPlus);

    //////////////////////////////////////////////////
    //Minus
    std::vector<particles> track_minus;
    track_minus.clear();
    for(int i = 0; i < total_particles.size(); i++)
        {
            particles tempParticle;
            tempParticle.arm          = total_particles[i].arm;
            tempParticle.charge       = total_particles[i].charge;
            tempParticle.energy       = total_particles[i].energy;
            tempParticle.eT           = total_particles[i].eT;
            tempParticle.eta          = total_particles[i].eta;
            tempParticle.phi          = total_particles[i].phi;
            tempParticle.phiDC        = total_particles[i].phiDC;
            tempParticle.ertTrigger   = total_particles[i].ertTrigger;
            tempParticle.id           = total_particles[i].id;
            tempParticle.jetPt        = total_particles[i].jetPt;

            if(tempParticle.charge != 0)
                {
                    float factor = 0.02;
                    if(total_particles[i].pT > 10)
                        {
                            factor = 0.02 + 0.01 * ((total_particles[i].pT - 10.0) / 10.0);
                        }
                    tempParticle.pT           = total_particles[i].pT - (factor * total_particles[i].pT);
                    tempParticle.mom          = tempParticle.pT / (total_particles[i].pT / total_particles[i].mom);
                    tempParticle.px           = tempParticle.pT * cos(total_particles[i].phi);
                    tempParticle.py           = tempParticle.pT * sin(total_particles[i].phi);
                    tempParticle.pz           = total_particles[i].pz;
                }
            else
                {
                    tempParticle.pT           = total_particles[i].pT;
                    tempParticle.mom          = total_particles[i].mom;
                    tempParticle.px           = total_particles[i].px;
                    tempParticle.py           = total_particles[i].py;
                    tempParticle.pz           = total_particles[i].pz;
                }
            track_minus.push_back(tempParticle);
        }
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, track_minus,
                 hPtRecoTrackMinus, hPtTrueMatchedTrackMinus, hPtRecoMatchedTrackMinus, hResponseMatrixTrackMinus);

    //////////////////////////////////////////////////
    //Cluster energy varied by +-3%
    //////////////////////////////////////////////////

    //////////////////////////////////////////////////
    //Plus
    std::vector<particles> cluster_plus;
    cluster_plus.clear();
    for(int i = 0; i < total_particles.size(); i++)
        {
            particles tempParticle;
            tempParticle.arm          = total_particles[i].arm;
            tempParticle.charge       = total_particles[i].charge;
            tempParticle.mom          = total_particles[i].mom;
            tempParticle.pT           = total_particles[i].pT;
            tempParticle.eta          = total_particles[i].eta;
            tempParticle.phi          = total_particles[i].phi;
            tempParticle.phiDC        = total_particles[i].phiDC;
            tempParticle.ertTrigger   = total_particles[i].ertTrigger;
            tempParticle.id           = total_particles[i].id;
            tempParticle.jetPt        = total_particles[i].jetPt;

            if(tempParticle.charge == 0)
                {
                    float factor = 0.03;
                    tempParticle.energy       = total_particles[i].energy + (factor * total_particles[i].energy);
                    tempParticle.eT           = tempParticle.energy * (total_particles[i].eT / total_particles[i].energy);
                    tempParticle.px           = tempParticle.eT * cos(total_particles[i].phi);
                    tempParticle.py           = tempParticle.eT * sin(total_particles[i].phi);
                    tempParticle.pz           = total_particles[i].pz;
                }
            else
                {
                    tempParticle.eT           = total_particles[i].eT;
                    tempParticle.energy       = total_particles[i].energy;
                    tempParticle.px           = total_particles[i].px;
                    tempParticle.py           = total_particles[i].py;
                    tempParticle.pz           = total_particles[i].pz;
                }
            cluster_plus.push_back(tempParticle);
        }
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, cluster_plus,
                 hPtRecoClusterPlus, hPtTrueMatchedClusterPlus, hPtRecoMatchedClusterPlus, hResponseMatrixClusterPlus);

    //////////////////////////////////////////////////
    //Minus
    std::vector<particles> cluster_minus;
    cluster_minus.clear();
    for(int i = 0; i < total_particles.size(); i++)
        {
            particles tempParticle;
            tempParticle.arm          = total_particles[i].arm;
            tempParticle.charge       = total_particles[i].charge;
            tempParticle.mom          = total_particles[i].mom;
            tempParticle.pT           = total_particles[i].pT;
            tempParticle.eta          = total_particles[i].eta;
            tempParticle.phi          = total_particles[i].phi;
            tempParticle.phiDC        = total_particles[i].phiDC;
            tempParticle.ertTrigger   = total_particles[i].ertTrigger;
            tempParticle.id           = total_particles[i].id;
            tempParticle.jetPt        = total_particles[i].jetPt;

            if(tempParticle.charge == 0)
                {
                    float factor = 0.03;
                    tempParticle.energy       = total_particles[i].energy - (factor * total_particles[i].energy);
                    tempParticle.eT           = tempParticle.energy * (total_particles[i].eT / total_particles[i].energy);
                    tempParticle.px           = tempParticle.eT * cos(total_particles[i].phi);
                    tempParticle.py           = tempParticle.eT * sin(total_particles[i].phi);
                    tempParticle.pz           = total_particles[i].pz;
                }
            else
                {
                    tempParticle.eT           = total_particles[i].eT;
                    tempParticle.energy       = total_particles[i].energy;
                    tempParticle.px           = total_particles[i].px;
                    tempParticle.py           = total_particles[i].py;
                    tempParticle.pz           = total_particles[i].pz;
                }
            cluster_minus.push_back(tempParticle);
        }
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, cluster_minus,
                 hPtRecoClusterMinus, hPtTrueMatchedClusterMinus, hPtRecoMatchedClusterMinus, hResponseMatrixClusterMinus);

    //////////////////////////////////////////////////
    //pT>2, energy>2
    //////////////////////////////////////////////////
    std::vector<particles> particles_tight;
    particles_tight.clear();
    for(int i = 0; i < total_particles.size(); i++)
        {
            int charge = total_particles[i].charge;
            float pT = total_particles[i].pT;
            float energy = total_particles[i].energy;

            bool passPt = (charge != 0) && (pT > 2.0);
            bool passEnergy = (charge == 0) && (energy > 2.0);

            if(passPt || passEnergy)
                {
                    particles_tight.push_back(total_particles[i]);
                }
        }
    getVariation(nc, minCf, maxCf, 0,
                 true_jets, particles_tight,
                 hPtRecoTrClTight, hPtTrueMatchedTrClTight, hPtRecoMatchedTrClTight, hResponseMatrixTrClTight);

    fillTrees(treesFill);

    if(nTotalEvents == 10)
        {
            float zVertexPythia = eventContainer->GetPrimaryVertexX();
            float zVertexPisaReco = eventContainer->GetPrimaryVertexY();
            float zVertexData = eventContainer->GetPrimaryVertexZ();

            cout << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Vertex info for 10th event: " << endl;
            cout << "PYTHIA: " << zVertexPythia << endl;
            cout << "PISA + RECO: " << zVertexPisaReco << endl;
            cout << "Data: " << zVertexData << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << endl;
        }

    return 0;
}

void JetAnalyzerSim::getPythiaTrueJet(std::vector<jets> true_jet_list,
                                      TH1F *hPtTrue[5],
                                      TH1F *hPtTrueEast[5],
                                      TH1F *hPtTrueWest[5])
{
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            int arm = true_jet_list[t].arm;
            float pT = true_jet_list[t].pT;

            hPtTrue[0]->Fill(pT);
            if(centralityBin == 1)
                {
                    hPtTrue[1]->Fill(pT);
                }
            if(centralityBin == 2)
                {
                    hPtTrue[2]->Fill(pT);
                }
            if(centralityBin == 3)
                {
                    hPtTrue[3]->Fill(pT);
                }
            if(centralityBin == 4)
                {
                    hPtTrue[4]->Fill(pT);
                }

            if(arm == 0)
                {
                    hPtTrueEast[0]->Fill(pT);
                    if(centralityBin == 1)
                        {
                            hPtTrueEast[1]->Fill(pT);
                        }
                    if(centralityBin == 2)
                        {
                            hPtTrueEast[2]->Fill(pT);
                        }
                    if(centralityBin == 3)
                        {
                            hPtTrueEast[3]->Fill(pT);
                        }
                    if(centralityBin == 4)
                        {
                            hPtTrueEast[4]->Fill(pT);
                        }
                }
            if(arm == 1)
                {
                    hPtTrueWest[0]->Fill(pT);
                    if(centralityBin == 1)
                        {
                            hPtTrueWest[1]->Fill(pT);
                        }
                    if(centralityBin == 2)
                        {
                            hPtTrueWest[2]->Fill(pT);
                        }
                    if(centralityBin == 3)
                        {
                            hPtTrueWest[3]->Fill(pT);
                        }
                    if(centralityBin == 4)
                        {
                            hPtTrueWest[4]->Fill(pT);
                        }
                }
        }
}


void JetAnalyzerSim::getRecoJet(std::vector<particles> total_particle_list,
                                std::vector<jets>& reco_jet_list,
                                TH1F *hPtReco[5])
{
    std::vector<particles> constituent_list;
    constituent_list.clear();
    jetAnalyzer->GetAntiKtCommon(total_particle_list, R, nc, minPt, minCf, maxCf, reco_jet_list, constituent_list);

    for(unsigned int s = 0; s < reco_jet_list.size(); s++)
        {
            float pT = reco_jet_list[s].pT;

            hPtReco[0]->Fill(pT);
            if(centralityBin == 1)
                {
                    hPtReco[1]->Fill(pT);
                }
            if(centralityBin == 2)
                {
                    hPtReco[2]->Fill(pT);
                }
            if(centralityBin == 3)
                {
                    hPtReco[3]->Fill(pT);
                }
            if(centralityBin == 4)
                {
                    hPtReco[4]->Fill(pT);
                }
        }
}


void JetAnalyzerSim::getMatchedJet(float fMinDeltaR,
                                   std::vector<jets> true_jet_list,
                                   std::vector<jets> reco_jet_list,
                                   std::vector<jets>& matched_jet_list,
                                   TH1F *hDistance[5],
                                   TH2F *hDistanceVsPtTrue[5],
                                   TH1F *hPtTrueMatched[5],
                                   TH1F *hPtRecoMatched[5],
                                   TH2F *hResponseMatrix[5],
                                   TH1F *hForJESR[5][22])
{
    std::vector<jetPair> jet_pair_list;
    jet_pair_list.clear();
    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            idTrue++;
            int tArm         = true_jet_list[t].arm;
            float tPt        = true_jet_list[t].pT;
            float tEta       = true_jet_list[t].eta;
            float tPhi       = true_jet_list[t].phi;
            float tNc        = true_jet_list[t].nc;
            float tCf        = true_jet_list[t].cf;
            float tNf        = true_jet_list[t].nf;
            float tDisc      = true_jet_list[t].disc;

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < reco_jet_list.size(); g++)
                {
                    idGood++;
                    float gArm       = reco_jet_list[g].arm;
                    float gPt        = reco_jet_list[g].pT;
                    float gEta       = reco_jet_list[g].eta;
                    float gPhi       = reco_jet_list[g].phi;
                    float gNc        = reco_jet_list[g].nc;
                    float gCf        = reco_jet_list[g].cf;
                    float gNf        = reco_jet_list[g].nf;
                    float gDisc      = reco_jet_list[g].disc;

                    if (tArm == gArm)
                        {
                            float deltaR = jetAnalyzer->dR(tEta, gEta, tPhi, gPhi);

                            jetPair tempPair;
                            tempPair.id.first    = idTrue;
                            tempPair.id.second   = idGood;
                            tempPair.pT.first    = tPt;
                            tempPair.pT.second   = gPt;
                            tempPair.eta.first   = tEta;
                            tempPair.eta.second  = gEta;
                            tempPair.phi.first   = tPhi;
                            tempPair.phi.second  = gPhi;
                            tempPair.nc.first    = tNc;
                            tempPair.nc.second   = gNc;
                            tempPair.cf.first    = tCf;
                            tempPair.cf.second   = gCf;
                            tempPair.nf.first    = tNf;
                            tempPair.nf.second   = gNf;
                            tempPair.disc.first  = tDisc;
                            tempPair.disc.second = gDisc;
                            tempPair.arm         = tArm;
                            tempPair.deltaR      = deltaR;
                            tempPair.centralityBin = centralityBin;

                            jet_pair_list.push_back(tempPair);
                        }
                }
        }

    //Sort pair by ascending order of deltaR
    std::sort(jet_pair_list.begin(), jet_pair_list.end(), sortPair());

    //Require 1 to 1 matching- and save as unique pair
    std::vector<jetPair> unique_pair_list;
    while (jet_pair_list.size())
        {
            unique_pair_list.push_back(jet_pair_list.front());
            jet_pair_list.erase(std::remove_if(jet_pair_list.begin(), jet_pair_list.end(), removePairId(jet_pair_list.front())),
                                jet_pair_list.end());
        }

    for(unsigned int u = 0; u < unique_pair_list.size(); u++)
        {
            float truePt = unique_pair_list[u].pT.first;
            float recoPt = unique_pair_list[u].pT.second;

            int arm      = unique_pair_list[u].arm;
            float deltaR = unique_pair_list[u].deltaR;

            hDistance[0]->Fill(deltaR);
            hDistanceVsPtTrue[0]->Fill(truePt, deltaR);
            if(centralityBin == 1)
                {
                    hDistance[1]->Fill(deltaR);
                    hDistanceVsPtTrue[1]->Fill(truePt, deltaR);
                }
            if(centralityBin == 2)
                {
                    hDistance[2]->Fill(deltaR);
                    hDistanceVsPtTrue[2]->Fill(truePt, deltaR);
                }
            if(centralityBin == 3)
                {
                    hDistance[3]->Fill(deltaR);
                    hDistanceVsPtTrue[3]->Fill(truePt, deltaR);
                }
            if(centralityBin == 4)
                {
                    hDistance[4]->Fill(deltaR);
                    hDistanceVsPtTrue[4]->Fill(truePt, deltaR);
                }

            if (deltaR < fMinDeltaR)
                {
                    hPtTrueMatched[0]->Fill(truePt);
                    hPtRecoMatched[0]->Fill(recoPt);
                    hResponseMatrix[0]->Fill(recoPt, truePt);
                    if(centralityBin == 1)
                        {
                            hPtTrueMatched[1]->Fill(truePt);
                            hPtRecoMatched[1]->Fill(recoPt);
                            hResponseMatrix[1]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 2)
                        {
                            hPtTrueMatched[2]->Fill(truePt);
                            hPtRecoMatched[2]->Fill(recoPt);
                            hResponseMatrix[2]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 3)
                        {
                            hPtTrueMatched[3]->Fill(truePt);
                            hPtRecoMatched[3]->Fill(recoPt);
                            hResponseMatrix[3]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 4)
                        {
                            hPtTrueMatched[4]->Fill(truePt);
                            hPtRecoMatched[4]->Fill(recoPt);
                            hResponseMatrix[4]->Fill(recoPt, truePt);
                        }

                    jets temp;
                    temp.arm           = unique_pair_list[u].arm;
                    temp.centralityBin = centralityBin;
                    temp.pT            = unique_pair_list[u].pT.second;
                    temp.eta           = unique_pair_list[u].eta.second;
                    temp.phi           = unique_pair_list[u].phi.second;
                    temp.nc            = unique_pair_list[u].nc.second;
                    temp.cf            = unique_pair_list[u].cf.second;
                    temp.nf            = unique_pair_list[u].nf.second;
                    temp.disc          = unique_pair_list[u].disc.second;

                    matched_jet_list.push_back(temp);

                    /////////////////////////////////////////////////////////////////////////////////////////////
                    //Plots for JES and JER
                    /////////////////////////////////////////////////////////////////////////////////////////////
                    float factor = recoPt / truePt;

                    fillJESR(hForJESR, 0, truePt, factor);
                    if(centralityBin == 1)
                        {
                            fillJESR(hForJESR, 1, truePt, factor);
                        }
                    if(centralityBin == 2)
                        {
                            fillJESR(hForJESR, 2, truePt, factor);
                        }
                    if(centralityBin == 3)
                        {
                            fillJESR(hForJESR, 3, truePt, factor);
                        }
                    if(centralityBin == 4)
                        {
                            fillJESR(hForJESR, 4, truePt, factor);
                        }
                }
        }
}

void JetAnalyzerSim::getVariation(float vNc, float vMinCf, float vMaxCf, int fiducialCut,
                                  std::vector<jets> true_jet_list,
                                  std::vector<particles> total_particle_list,
                                  TH1F *hPtReco[5],
                                  TH1F *hPtTrueMatched[5],
                                  TH1F *hPtRecoMatched[5],
                                  TH2F *hResponseMatrix[5])
{
    std::vector<particles> constituent_list;
    constituent_list.clear();
    std::vector<jets> reco_jet_list;
    reco_jet_list.clear();


    jetAnalyzer->GetAntiKtCommon(total_particle_list, R, vNc, minPt, vMinCf, vMaxCf, reco_jet_list, constituent_list, fiducialCut);

    for(unsigned int g = 0; g < reco_jet_list.size(); g++)
        {
            float recoPt = reco_jet_list[g].pT;

            hPtReco[0]->Fill(recoPt);
            if(centralityBin == 1)
                {
                    hPtReco[1]->Fill(recoPt);
                }
            if(centralityBin == 2)
                {
                    hPtReco[2]->Fill(recoPt);
                }
            if(centralityBin == 3)
                {
                    hPtReco[3]->Fill(recoPt);
                }
            if(centralityBin == 4)
                {
                    hPtReco[4]->Fill(recoPt);
                }
        }

    std::vector<jetPair> jet_pair_list;
    jet_pair_list.clear();
    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            idTrue++;
            int tArm         = true_jet_list[t].arm;
            float tPt        = true_jet_list[t].pT;
            float tEta       = true_jet_list[t].eta;
            float tPhi       = true_jet_list[t].phi;
            float tNc        = true_jet_list[t].nc;
            float tCf        = true_jet_list[t].cf;
            float tNf        = true_jet_list[t].nf;
            float tDisc      = true_jet_list[t].disc;

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < reco_jet_list.size(); g++)
                {
                    idGood++;
                    float gArm       = reco_jet_list[g].arm;
                    float gPt        = reco_jet_list[g].pT;
                    float gEta       = reco_jet_list[g].eta;
                    float gPhi       = reco_jet_list[g].phi;
                    float gNc        = reco_jet_list[g].nc;
                    float gCf        = reco_jet_list[g].cf;
                    float gNf        = reco_jet_list[g].nf;
                    float gDisc      = reco_jet_list[g].disc;

                    if (tArm == gArm)
                        {
                            float deltaR = jetAnalyzer->dR(tEta, gEta, tPhi, gPhi);

                            jetPair tempPair;
                            tempPair.id.first    = idTrue;
                            tempPair.id.second   = idGood;
                            tempPair.pT.first    = tPt;
                            tempPair.pT.second   = gPt;
                            tempPair.eta.first   = tEta;
                            tempPair.eta.second  = gEta;
                            tempPair.phi.first   = tPhi;
                            tempPair.phi.second  = gPhi;
                            tempPair.nc.first    = tNc;
                            tempPair.nc.second   = gNc;
                            tempPair.cf.first    = tCf;
                            tempPair.cf.second   = gCf;
                            tempPair.nf.first    = tNf;
                            tempPair.nf.second   = gNf;
                            tempPair.disc.first  = tDisc;
                            tempPair.disc.second = gDisc;
                            tempPair.arm         = tArm;
                            tempPair.deltaR      = deltaR;
                            tempPair.centralityBin = centralityBin;

                            jet_pair_list.push_back(tempPair);
                        }
                }
        }

    //Sort pair by ascending order of deltaR
    std::sort(jet_pair_list.begin(), jet_pair_list.end(), sortPair());

    //Require 1 to 1 matching- and save as unique pair
    std::vector<jetPair> unique_pair_list;
    while (jet_pair_list.size())
        {
            unique_pair_list.push_back(jet_pair_list.front());
            jet_pair_list.erase(std::remove_if(jet_pair_list.begin(), jet_pair_list.end(), removePairId(jet_pair_list.front())),
                                jet_pair_list.end());
        }

    for(unsigned int u = 0; u < unique_pair_list.size(); u++)
        {
            float truePt = unique_pair_list[u].pT.first;
            float recoPt = unique_pair_list[u].pT.second;

            int arm      = unique_pair_list[u].arm;
            float deltaR = unique_pair_list[u].deltaR;

            if (deltaR < minDeltaR)
                {
                    hPtTrueMatched[0]->Fill(truePt);
                    hPtRecoMatched[0]->Fill(recoPt);
                    hResponseMatrix[0]->Fill(recoPt, truePt);
                    if(centralityBin == 1)
                        {
                            hPtTrueMatched[1]->Fill(truePt);
                            hPtRecoMatched[1]->Fill(recoPt);
                            hResponseMatrix[1]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 2)
                        {
                            hPtTrueMatched[2]->Fill(truePt);
                            hPtRecoMatched[2]->Fill(recoPt);
                            hResponseMatrix[2]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 3)
                        {
                            hPtTrueMatched[3]->Fill(truePt);
                            hPtRecoMatched[3]->Fill(recoPt);
                            hResponseMatrix[3]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 4)
                        {
                            hPtTrueMatched[4]->Fill(truePt);
                            hPtRecoMatched[4]->Fill(recoPt);
                            hResponseMatrix[4]->Fill(recoPt, truePt);
                        }
                }
        }
}


void JetAnalyzerSim::fillTrees(bool what)
{
    if(what)
        {
            tTrueJets->Fill();

            tTotalParticles->Fill();

            tRecoJets->Fill();
            tMatchedJets->Fill();
        }
}


int JetAnalyzerSim::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetAnalyzerSim::End(PHCompositeNode *topNode)
{
    //Bins 1 to 5: Total Events
    hEvents->SetBinContent(1, nTotalEvents);
    hEvents->SetBinContent(2, nTotalEventsCentrality1);
    hEvents->SetBinContent(3, nTotalEventsCentrality2);
    hEvents->SetBinContent(4, nTotalEventsCentrality3);
    hEvents->SetBinContent(5, nTotalEventsCentrality4);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nTotalEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}

























































