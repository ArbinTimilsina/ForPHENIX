//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//PHPythia tools
#include <PHPyTrigger.h>
#include <TMCParticle.h>
#include <PHPythiaHeader.h>
#include <PHPyCommon.h>
#include <PHPythiaContainer.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TNtuple.h>


//My source file
#include "Run13JetSim.h"

using namespace std;
using namespace findNode;
using namespace fastjet;
//================================ Constructor ================================
//Here we can initiate some variables
Run13JetSim::Run13JetSim(float iminPt)
    : SubsysReco("Run13JetSim")
{
    minPt = iminPt;

    nTotalEvents = 0;

    treesFill = false;

    return;
}

int Run13JetSim::Init(PHCompositeNode *topNode)
{

    //To access various nodes
    se = Fun4AllServer::instance();
    rc = recoConsts::instance();

    //For PYTHIA trigger
    run13JetTriggerPythia = new Run13JetTriggerPythia(minPt, 1);

    //For tracks and clusters
    run13Jet = new Run13Jet("dummyFile.root");
    run13Jet->Init(topNode);

    //Output file name
    outfile = new TFile("Run13JetSim.root", "RECREATE");

    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    hVertexPythia = new TH1F("hVertexPythia", "Vertex Distribution for PYTHIA", 200, -50, 50);
    hVertexReco = new TH1F("hVertexReco", "Vertex Distribution for Reco", 200, -50, 50);
    hDistance = new TH1F("hDistance", "#DeltaR between True Jet and Reco Jet", 400, 0.0, 0.4);

    hPtTrueJet = new TH1F("hPtTrueJet", "p_{T} of True Jet",  NPTBINS_TRUE, PTBINS_TRUE);
    hPtTrueMatchedJet = new TH1F("hPtTrueMatchedJet", "p_{T} of True Jet Matched",  NPTBINS_TRUE, PTBINS_TRUE);

    hPtRecoJet = new TH1F("hPtRecoJet", "p_{T} of Reco Jet",  NPTBINS_RECO, PTBINS_RECO);
    hPtRecoMatchedJet = new TH1F("hPtRecoMatchedJet", "p_{T} of Reco Jet Matched",  NPTBINS_RECO, PTBINS_RECO);

    hResponseMatrix = new TH2F("hResponseMatrix", "p_{T, true} vs p_{T, reco} for Matched Jets",
                               NPTBINS_RECO, PTBINS_RECO, NPTBINS_TRUE, PTBINS_TRUE);


    hEMCalTOF = new TH1F("hEMCalTOF", "EMCal TOF", 400, -100, 100);
    for (int ias = 0; ias < NARMSECT; ias++)
        {
            int ny, nz;
            if (IsPbGl(ias))
                {
                    ny = YPOS_PBGL;
                    nz = ZPOS_PBGL;
                }
            else
                {
                    ny = YPOS_PBSC;
                    nz = ZPOS_PBSC;
                }
            hSectorHits[ias] = new TH2F(Form("hSecorHits_%i", ias), Form("Sector %i", ias), nz, 0, nz, ny, 0, ny);
        }
    hClustersEta = new TH1F("hClustersEta", "Clusters going into Jet reconstruction", 160, -0.4, 0.4);
    hClustersPhi = new TH1F("hClustersPhi", "Clusters going into Jet reconstruction", 200, -1.0, 4.0);

    hModifiedQuality_NE = new TH2F("hModifiedQuality_NE", "Modified Quality Cut, NE",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_SE = new TH2F("hModifiedQuality_SE", "Modified Quality Cut, SE",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_NW = new TH2F("hModifiedQuality_NW", "Modified Quality Cut, NW",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_SW = new TH2F("hModifiedQuality_SW", "Modified Quality Cut, SW",  400, 0, 80, 120, -0.6, 0.6);

    hTracksZed = new TH1F("hTracksZed", "Tracks going into Jet reconstruction", 400, -100.0, 100.0);
    hTracksPhi = new TH1F("hTracksPhi",	"Tracks going into Jet reconstruction", 200, -1.0, 4.0);
    //*******************************************************************************************************

    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    if(treesFill)
        {
            tPythiaParticle = new TTree("tPythiaParticle", "Final particles of PYTHIA");
            tPythiaParticle->Branch("tPythiaParticle", &pythia_particle);
            tPythiaTrueJet = new TTree("tPythiaTrueJet", "True Jet of Pythia");
            tPythiaTrueJet->Branch("tPythiaTrueJet", &pythia_true_jet);

            tRecoTrack = new TTree("tRecoTrack", "Simulated good tracks");
            tRecoTrack->Branch("tRecoTrack", &charged_particle);
            tRecoCluster = new TTree("tRecoCluster", "Simulated good clusters");
            tRecoCluster->Branch("tRecoCluster", &neutral_particle);

            tRecoJet = new TTree("tRecoJet", "Reco jet");
            tRecoJet->Branch("tRecoJet", &reco_jet);
            tMatchedJet = new TTree("tMatchedJet", "Matched jet");
            tMatchedJet->Branch("tMatchedJet", &matched_jet);
        }

    return EVENT_OK;
}

int Run13JetSim::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int Run13JetSim::ResetEvent(PHCompositeNode *topNode)
{
    pythia_particle.clear();
    pythia_true_jet.clear();

    charged_particle.clear();
    neutral_particle.clear();

    reco_jet.clear();
    matched_jet.clear();

    return EVENT_OK;
}


int Run13JetSim::process_event(PHCompositeNode *topNode)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Get True Jet
    ///////////////////////////////////////////////////////////////////////////////////////////////
    getPythiaTrueJet(topNode, pythia_particle, pythia_true_jet, hVertexPythia, hPtTrueJet);

    //After vertex feeding on the PYTHIA particles, true jet might not reconstruct
    if(pythia_true_jet.size() == 0)
        {
            return ABORTEVENT;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Get Simulated tracks and clusters
    ///////////////////////////////////////////////////////////////////////////////////////////////
    PHCompositeNode *pisaRecoNode = se->topNode(rc->get_CharFlag("PISARECO_TOPNODE"));
    PHGlobal *phglobal = getClass<PHGlobal>(pisaRecoNode, "PHGlobal");
    if (!phglobal)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    float zvertex = phglobal->getBbcZVertex();
    if(nTotalEvents == 1)
        {
            cout << "Vertex of reco: " << zvertex << endl;
            cout << "*********************************" << endl << endl;
        }
    hVertexReco->Fill(zvertex);

    float bbcT0 = phglobal->getBbcTimeZero();
    run13Jet->GetParticles(pisaRecoNode, charged_particle, neutral_particle, zvertex, bbcT0, false);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //Make anti-kt jet here
    ////////////////////////////////////////////////////////////////////////////////////////////////
    float R = 0.3;
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    const int indexTotal = charged_particle.size() + neutral_particle.size();
    unsigned int indexCharged = charged_particle.size();

    float particlePt[indexTotal];
    fill(particlePt, particlePt + indexTotal / sizeof(float), -999.9);

    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    int index = 0;
    for (unsigned int h = 0; h < charged_particle.size(); h++)
        {
            fastjet::PseudoJet pseudoCharged(charged_particle[h].px,
                                             charged_particle[h].py,
                                             charged_particle[h].pz,
                                             charged_particle[h].mom);
            pseudoCharged.set_user_index(index);
            particlePt[index] = charged_particle[h].pT;
            jetParticles_all.push_back(pseudoCharged);
            index++;
        }

    for (unsigned int n = 0; n < neutral_particle.size(); n++)
        {
            fastjet::PseudoJet pseudoNeutral(neutral_particle[n].px,
                                             neutral_particle[n].py,
                                             neutral_particle[n].pz,
                                             neutral_particle[n].energy);
            pseudoNeutral.set_user_index(index);
            particlePt[index] = neutral_particle[n].pT;
            jetParticles_all.push_back(pseudoNeutral);
            index++;
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();

    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float chargedPt    = 0.0;

            recoJets tempJet;
            tempJet.pT = aFastJet.perp();
            tempJet.eta = aFastJet.pseudorapidity();
            tempJet.phi =  run13JetTriggerPythia->phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();
            for (unsigned int iconst = 0; iconst < nconst; iconst++)
                {
                    unsigned int indx = constituents[iconst].user_index();

                    if (indx < indexCharged)//Charged particles
                        {
                            chargedPt += particlePt[indx];
                        }
                }
            tempJet.nc = (float)nconst;
            tempJet.cf = chargedPt / tempJet.pT;

            bool passJetLevelCuts = (tempJet.pT > minPt) && (tempJet.nc >= 3.0) && (tempJet.cf > 0.2) && (tempJet.cf < 0.7);
            if(passJetLevelCuts)
                {
                    reco_jet.push_back(tempJet);
                    hPtRecoJet->Fill(tempJet.pT);
                }
        }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Do matching here
    ///////////////////////////////////////////////////////////////////////////////////////////////
    getMatchedJet(pythia_true_jet, reco_jet, matched_jet, hDistance, hPtTrueMatchedJet, hPtRecoMatchedJet, hResponseMatrix);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Fill extra histograms
    ///////////////////////////////////////////////////////////////////////////////////////////////
    for (unsigned int n = 0; n < neutral_particle.size(); n++)
        {
            hEMCalTOF->Fill(neutral_particle[n].tof);

            int armsect = neutral_particle[n].armsect;
            int yTowerPos = neutral_particle[n].yTowerPos;
            int zTowerPos = neutral_particle[n].zTowerPos;
            hSectorHits[armsect]->Fill(zTowerPos, yTowerPos);

            hClustersEta->Fill(neutral_particle[n].eta);
            hClustersPhi->Fill(neutral_particle[n].phi);
        }

    for (unsigned int c = 0; c < charged_particle.size(); c++)
        {
            hTracksZed->Fill(charged_particle[c].zedDC);
            hTracksPhi->Fill(charged_particle[c].phiDC);

            float pT = charged_particle[c].pT;
            float zedDC = charged_particle[c].zedDC;
            int arm = charged_particle[c].arm;

            float board = charged_particle[c].board;
            float alpha = charged_particle[c].alpha;

            if (pT > 0.2 && pT < 25)
                {
                    if (zedDC > 0 && arm == 0)
                        {
                            hModifiedQuality_NE->Fill(board, alpha);
                        }
                    if (zedDC > 0 && arm == 1)
                        {
                            hModifiedQuality_NW->Fill(board, alpha);
                        }
                    if (zedDC < 0 && arm == 0)
                        {
                            hModifiedQuality_SE->Fill(board, alpha);
                        }
                    if (zedDC < 0 && arm == 1)
                        {
                            hModifiedQuality_SW->Fill(board, alpha);
                        }
                }
        }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    fillTrees(treesFill);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    nTotalEvents++;

    return 0;
}

void Run13JetSim::getPythiaTrueJet(PHCompositeNode *topNode,
                                   std::vector<pythiaParticles>& particle_list,
                                   std::vector<pythiaJets>& true_jet_list,
                                   TH1F *hVertex, TH1F *hPtTrue)
{
    PHCompositeNode *pythiaNode = se->topNode(rc->get_CharFlag("PHPYTHIA_TOPNODE"));
    PHPythiaHeader *phpythiaheader = findNode::getClass<PHPythiaHeader>(pythiaNode, "PHPythiaHeader");
    if (!phpythiaheader)
        {
            cout << "No PHPythiaHeader! No sense continuing" << endl;
            exit(1);
        }

    float zvertex = phpythiaheader->GetPrimaryVertexZ();
    zvertex = zvertex * MM2CM;
    hVertex->Fill(zvertex);

    if(nTotalEvents == 1)
        {
            cout << endl;
            cout << "*********************************" << endl;
            cout << "Vertex matching info for event 1: " << endl;
            cout << "*********************************" << endl;
            cout << "Vertex of Pythia: " << zvertex << endl;
        }

    run13JetTriggerPythia->GetPythia(pythiaNode, particle_list, true_jet_list);
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            hPtTrue->Fill(true_jet_list[t].pT);
        }
}

void Run13JetSim::getMatchedJet(std::vector<pythiaJets> true_jet_list, std::vector<recoJets> reco_jet_list, std::vector<recoJets>& matched_jet_list,
                                TH1F *hDis, TH1F *hPtTrueM, TH1F *hPtRecoM, TH2F *hResMatrix)
{
    std::vector<jetPair> jet_pair_list;
    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            idTrue++;
            float tPt = true_jet_list[t].pT;
            float tEta = true_jet_list[t].eta;
            float tPhi = true_jet_list[t].phi;
            float tNc = -9999.9;
            float tCf = -9999.9;

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < reco_jet_list.size(); g++)
                {
                    idGood++;
                    float gPt =  reco_jet_list[g].pT;
                    float gEta = reco_jet_list[g].eta;
                    float gPhi = reco_jet_list[g].phi;
                    float gNc = reco_jet_list[g].nc;
                    float gCf = reco_jet_list[g].cf;

                    if (getArm(true_jet_list[t].phi) == getArm(reco_jet_list[g].phi))
                        {
                            float deltaR = dR(tEta, gEta, tPhi, gPhi);

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
                            tempPair.deltaR      = deltaR;

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
            float deltaR = unique_pair_list[u].deltaR;

            hDis->Fill(deltaR);

            if (deltaR < 0.3)
                {
                    hPtTrueM->Fill(truePt);
                    hPtRecoM->Fill(recoPt);
                    hResMatrix->Fill(recoPt, truePt);

                    recoJets temp;
                    temp.pT = unique_pair_list[u].pT.second;
                    temp.eta = unique_pair_list[u].eta.second;
                    temp.phi = unique_pair_list[u].phi.second;
                    temp.nc = unique_pair_list[u].nc.second;
                    temp.cf = unique_pair_list[u].cf.second;

                    matched_jet_list.push_back(temp);
                }
        }
}

void Run13JetSim::fillTrees(bool write)
{
    if(write)
        {
            tPythiaParticle->Fill();
            tPythiaTrueJet->Fill();

            tRecoTrack->Fill();
            tRecoCluster->Fill();

            tRecoJet->Fill();
            tMatchedJet->Fill();
        }
}

int Run13JetSim::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int Run13JetSim::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nTotalEvents);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nTotalEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}

















