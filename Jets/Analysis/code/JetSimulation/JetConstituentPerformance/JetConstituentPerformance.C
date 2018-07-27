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
#include <PHPythiaHeader.h>
#include <PHPyCommon.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TNtuple.h>

//For vertex during sim reco
#include <VtxOut.h>

//My source file
#include "JetConstituentPerformance.h"

//For acceptance
#include <TrackQualityPP.h>

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
JetConstituentPerformance::JetConstituentPerformance(const float ir, const float iminPt,
						     const float incPythia, const float iminCfPythia, const float imaxCfPythia,
						     const float iminChargedDeltaR, const float iminNeutralDeltaR,
						     const bool iperfectEMCal, const bool iphenixParticle)
    : SubsysReco("JetConstituentPerformance")
{
    R                 = ir;
    minPt             = iminPt;

    ncPythia          = incPythia;
    minCfPythia       = iminCfPythia;
    maxCfPythia       = imaxCfPythia;

    minChargedDeltaR  = iminChargedDeltaR;
    minNeutralDeltaR  = iminNeutralDeltaR;

    perfectEMCal      = iperfectEMCal;
    phenixParticle    = iphenixParticle;

    nTotalEvents      = 0;

    nTotalPyParticles = 0;

    treesFill         = false;

    return;
}

int JetConstituentPerformance::Init(PHCompositeNode *topNode)
{

    //To access various nodes
    se = Fun4AllServer::instance();
    rc = recoConsts::instance();

    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(false);
    jetAnalyzer->SetCuAu(false);
    jetAnalyzer->MyInit();

    //For PYTHIA trigger
    jetTriggerPythia = new JetTriggerPythia(R, ncPythia, minPt, minCfPythia, maxCfPythia, 1, phenixParticle);

    //Output file name
    outfile = new TFile("JetConstituentPerformance.root", "RECREATE");

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    //*******************************************************************************************************
    hPtTrueJet = new TH1F("hPtTrueJet", "p_{T} of True Jets of PYTHIA",
                          NPTBINS_FINAL, PTBINS_FINAL);

    hPtPion = new TH1F("hPtPion", "p_{T} of Pions in True Jets of PYTHIA",
                       NPTBINS_FINAL, PTBINS_FINAL);

    hPtPhoton = new TH1F("hPtPhoton", "p_{T} of Photons in True Jets of PYTHIA",
                         NPTBINS_FINAL, PTBINS_FINAL);

    hPtKaon = new TH1F("hPtKaon", "p_{T} of Kaons in True Jets of PYTHIA",
                       NPTBINS_FINAL, PTBINS_FINAL);

    hPtProton = new TH1F("hPtProton", "p_{T} of Protons in True Jets of PYTHIA",
                         NPTBINS_FINAL, PTBINS_FINAL);

    hPtElectron = new TH1F("hPtElectron", "p_{T} of Electrons in True Jets of PYTHIA",
                           NPTBINS_FINAL, PTBINS_FINAL);

    hPtNeutron = new TH1F("hPtNeutron", "p_{T} of Neutrons in True Jets of PYTHIA",
                          NPTBINS_FINAL, PTBINS_FINAL);

    hPtK0L = new TH1F("hPtK0L", "p_{T} of K0L in True Jets of PYTHIA",
                      NPTBINS_FINAL, PTBINS_FINAL);

    //*******************************************************************************************************
    for (unsigned int c = 0; c < 2; c++)
        {
            hPEtTrueParticle[c] = new TH1F(Form("hPEtTrueParticle_%u", c),
                                           "p_{T} or E_{T} of Constituents of True Jet",
                                           20 * 2, 0.5, 20.5);

            hPEtRecoParticle[c] = new TH1F(Form("hPEtRecoParticle_%u", c),
                                           "p_{T} or E_{T} of Good Reconstructed Tracks and Clusters",
                                           20 * 2, 0.5, 20.5);

            hPEtForEfficiency[c] = new TH1F(Form("hPEtForEfficiency_%u", c),
                                            "p_{T} or E_{T} for Efficiency",
                                            20 * 2, 0.5, 20.5);

            hDistance[c] = new TH1F(Form("hDistance_%u", c), "#DeltaR between particles, and reco tracks and clusters",
                                    800, 0.0, 2.0);

            hResponseMatrix[c] = new TH2F(Form("hResponseMatrix_%u", c),
                                          "Response Matrix",
                                          20 * 2, 0.5, 20.5, 20 * 2, 0.5, 20.5);
        }

    //*******************************************************************************************************
    for (unsigned int c = 0; c < 4; c++)
        {
            hTracksDistance[c] = new TH1F(Form("hTracksDistance_%u", c), "#DeltaR between particles, and reco tracks with different cuts",
                                          800, 0.0, 2.0);

            hTracksEfficiency[c] = new TH1F(Form("hTracksEfficiency_%u", c),
                                            "p_{T} for Efficiency with various track cuts",
                                            20 * 2, 0.5, 20.5);
        }
    //*******************************************************************************************************

    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    if(treesFill)
        {
            tPythiaParticle = new TTree("tPythiaParticle", "Final particles of PYTHIA in Run 12 PHENIX acceptance");
            tPythiaParticle->Branch("tPythiaParticle", &pythia_particle_list);

            tPythiaTrueJet = new TTree("tPythiaTrueJet", "True Jet of Pythia");
            tPythiaTrueJet->Branch("tPythiaTrueJet", &pythia_true_jet_list);

            tPythiaTrueJetConstituent = new TTree("tPythiaTrueJetConstituent", "Constituents of PYTHIA True Jet");
            tPythiaTrueJetConstituent->Branch("tPythiaTrueJetConstituent", &pythia_true_jet_constituent_list);

            //*******************************************************************************************************
            tPisaRecoAllTrack = new TTree("tPisaRecoAllTrack", "Simulated All Tracks");
            tPisaRecoAllTrack->Branch("tPisaRecoAllTrack", &pisa_reco_all_track_list);

            tPisaRecoAllCluster = new TTree("tPisaRecoAllCluster", "Simulated All Clusters");
            tPisaRecoAllCluster->Branch("tPisaRecoAllCluster", &pisa_reco_all_cluster_list);

            tPisaRecoGoodTrack = new TTree("tPisaRecoGoodTrack", "Simulated Tracks with all cuts");
            tPisaRecoGoodTrack->Branch("tPisaRecoGoodTrack", &pisa_reco_good_track_list);

            tPisaRecoGoodCluster = new TTree("tPisaRecoGoodCluster", "Simulated Clusters with all cuts");
            tPisaRecoGoodCluster->Branch("tPisaRecoGoodCluster", &pisa_reco_good_cluster_list);

            tPisaRecoGoodParticle = new TTree("tPisaRecoGoodParticle", "Simulated Tracks and Clusters with all cuts");
            tPisaRecoGoodParticle->Branch("tPisaRecoGoodParticle", &pisa_reco_good_particle_list);
        }

    return EVENT_OK;
}

int JetConstituentPerformance::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int JetConstituentPerformance::ResetEvent(PHCompositeNode *topNode)
{
    pythia_particle_list.clear();
    pythia_true_jet_list.clear();
    pythia_true_jet_constituent_list.clear();

    pisa_reco_all_track_list.clear();
    pisa_reco_all_cluster_list.clear();

    pisa_reco_good_track_list.clear();
    pisa_reco_good_cluster_list.clear();
    pisa_reco_good_particle_list.clear();

    return EVENT_OK;
}


int JetConstituentPerformance::process_event(PHCompositeNode *topNode)
{
    //Get True Jet
    getPythiaTrueJet(topNode,
                     pythia_particle_list, pythia_true_jet_list, pythia_true_jet_constituent_list);

    nTotalPyParticles += pythia_particle_list.size();

    //Get stuff
    getEfficiency(topNode,
                  pythia_true_jet_constituent_list,
                  pisa_reco_all_track_list, pisa_reco_all_cluster_list,
                  pisa_reco_good_track_list, pisa_reco_good_cluster_list, pisa_reco_good_particle_list,
                  hPEtTrueParticle, hDistance, hPEtRecoParticle, hPEtForEfficiency, hResponseMatrix, perfectEMCal);

    //Get Track Reconstruction Efficiency
    getTrackCutsEfficiency(topNode,
                           pythia_true_jet_constituent_list,
                           hTracksDistance, hTracksEfficiency);

    ///////////////////////////////////////////////////////////////////////////////////////////////

    fillHistograms();

    fillTrees(treesFill);

    nTotalEvents++;

    return 0;
}


void JetConstituentPerformance::getPythiaTrueJet(PHCompositeNode *topNode,
						 std::vector<particles>& particle_list,
						 std::vector<jets>& true_jet_list,
						 std::vector<particles>& true_jet_constituent_list)
{
    PHCompositeNode *pythiaNode = se->topNode(rc->get_CharFlag("PHPYTHIA_TOPNODE"));

    if(nTotalEvents == 0)
        {
            jetTriggerPythia->GetPythia(pythiaNode, particle_list, true_jet_list, true_jet_constituent_list, true);

            cout << "For first triggered event: " << endl;
            cout << "Particle's:" << setw(15) << "charge" << setw(15) << "energy" << setw(15) << "mom" << setw(15)
                 << "pT" << setw(15) << "px" << setw(15) << "py" << setw(15) << "pz" << setw(15) << "eta" << setw(15) << "phi" << endl;
            for(int p = 0; p < particle_list.size(); p++)
                {
                    cout << "Particle #: " << p + 1 << setw(15) << particle_list[p].charge << setw(15) << particle_list[p].energy << setw(15)
                         << particle_list[p].mom << setw(15) << particle_list[p].pT << setw(15) << particle_list[p].px << setw(15)
                         << particle_list[p].py << setw(15) << particle_list[p].pz << setw(15) << particle_list[p].eta << setw(15)
                         << particle_list[p].phi << endl;
                }
            cout << "***********************************************************************" << endl << endl;
        }
    else
        {
            jetTriggerPythia->GetPythia(pythiaNode, particle_list, true_jet_list, true_jet_constituent_list, false);
        }
}


void JetConstituentPerformance::getEfficiency(PHCompositeNode *topNode,
					      std::vector<particles> true_jet_constituent_list,
					      std::vector<tracks>& all_track_list,
					      std::vector<clusters>& all_cluster_list,
					      std::vector<tracks>& good_track_list,
					      std::vector<clusters>& good_cluster_list,
					      std::vector<particles>& good_particle_list,
					      TH1F *hPEtTrue[2],
					      TH1F *hDeltaR[2],
					      TH1F *hPEtReco[2],
					      TH1F *hPEtEfficiency[2],
					      TH2F *hMatrix[2],
					      bool perfectEMC)
{
    PHCompositeNode *pisaRecoNode = se->topNode(rc->get_CharFlag("PISARECO_TOPNODE"));

    //For ZVertex used in reconstruction
    //https://www.phenix.bnl.gov/viewvc/viewvc.cgi/phenix/offline/packages/vtx/VtxOut.h?revision=1.25&view=co
    VtxOut *vertexOut = getClass<VtxOut>(pisaRecoNode, "VtxOut");
    if (!vertexOut)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    float zvertex = vertexOut->get_ZVertex();

    //Get Tracks
    jetAnalyzer->GetTracks(pisaRecoNode, all_track_list);

    //Get Clusters
    jetAnalyzer->GetClusters(pisaRecoNode, all_cluster_list, zvertex, perfectEMC);

    jetAnalyzer->GetParticles(all_track_list, all_cluster_list, good_track_list, good_cluster_list, good_particle_list);

    std::vector<particlePair> charged_particle_pair_list;
    charged_particle_pair_list.clear();
    std::vector<particlePair> neutral_particle_pair_list;
    neutral_particle_pair_list.clear();

    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_constituent_list.size(); t++)
        {
            idTrue++;
            float tPt        = true_jet_constituent_list[t].pT;
            float tEt        = true_jet_constituent_list[t].eT;
            float tEta       = true_jet_constituent_list[t].eta;
            float tPhi       = true_jet_constituent_list[t].phi;

            if(true_jet_constituent_list[t].charge == 0)
                {
                    hPEtTrue[0]->Fill(tEt);
                }
            else
                {
                    hPEtTrue[1]->Fill(tPt);
                }

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < good_particle_list.size(); g++)
                {
                    idGood++;
                    float gPt        = good_particle_list[g].pT;
                    float gEt        = good_particle_list[g].eT;
                    float gEta       = good_particle_list[g].eta;
                    float gPhi       = good_particle_list[g].phi;

                    if((true_jet_constituent_list[t].arm == good_particle_list[g].arm) &&
		       (true_jet_constituent_list[t].charge == good_particle_list[g].charge))
                        {
                            float deltaR = jetAnalyzer->dR(tEta, gEta, tPhi, gPhi);

                            particlePair tempPair;
                            tempPair.id.first  = idTrue;
                            tempPair.id.second = idGood;
                            tempPair.pT.first  = tPt;
                            tempPair.pT.second = gPt;
                            tempPair.eT.first  = tEt;
                            tempPair.eT.second = gEt;
                            tempPair.deltaR    = deltaR;

                            if(true_jet_constituent_list[t].charge == 0)
                                {
                                    neutral_particle_pair_list.push_back(tempPair);

                                }
                            else
                                {
                                    charged_particle_pair_list.push_back(tempPair);
                                }
                        }
                }
        }

    //Sort pair by ascending order of deltaR
    std::sort(charged_particle_pair_list.begin(), charged_particle_pair_list.end(), sortPair());
    std::sort(neutral_particle_pair_list.begin(), neutral_particle_pair_list.end(), sortPair());

    //Require 1 to 1 matching- and save as unique pair
    std::vector<particlePair> unique_charged_particle_pair_list;
    unique_charged_particle_pair_list.clear();
    while (charged_particle_pair_list.size())
        {
            unique_charged_particle_pair_list.push_back(charged_particle_pair_list.front());
            charged_particle_pair_list.erase(std::remove_if(charged_particle_pair_list.begin(), charged_particle_pair_list.end(),
							    removePairId(charged_particle_pair_list.front())), charged_particle_pair_list.end());
        }

    std::vector<particlePair> unique_neutral_particle_pair_list;
    unique_neutral_particle_pair_list.clear();
    while (neutral_particle_pair_list.size())
        {
            unique_neutral_particle_pair_list.push_back(neutral_particle_pair_list.front());
            neutral_particle_pair_list.erase(std::remove_if(neutral_particle_pair_list.begin(), neutral_particle_pair_list.end(),
							    removePairId(neutral_particle_pair_list.front())), neutral_particle_pair_list.end());
        }


    for(unsigned int u = 0; u < unique_neutral_particle_pair_list.size(); u++)
        {
            float trueEt = unique_neutral_particle_pair_list[u].eT.first;
            float recoEt = unique_neutral_particle_pair_list[u].eT.second;
            float deltaR = unique_neutral_particle_pair_list[u].deltaR;
            hDeltaR[0]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hPEtReco[0]->Fill(recoEt);
                    hPEtEfficiency[0]->Fill(trueEt);
                    hMatrix[0]->Fill(trueEt, recoEt);
                }
        }

    for(unsigned int u = 0; u < unique_charged_particle_pair_list.size(); u++)
        {
            float truePt = unique_charged_particle_pair_list[u].pT.first;
            float recoPt = unique_charged_particle_pair_list[u].pT.second;
            float deltaR = unique_charged_particle_pair_list[u].deltaR;
            hDeltaR[1]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hPEtReco[1]->Fill(recoPt);
                    hPEtEfficiency[1]->Fill(truePt);
                    hMatrix[1]->Fill(truePt, recoPt);
                }
        }
}

void JetConstituentPerformance::getTrackCutsEfficiency(PHCompositeNode *topNode,
						       std::vector<particles> true_jet_constituent_list,
						       TH1F *hTDistance[4],
						       TH1F *hTEfficiency[4])
{
    PHCompositeNode *pisaRecoNode = se->topNode(rc->get_CharFlag("PISARECO_TOPNODE"));

    std::vector<tracks> all_track_list;
    all_track_list.clear();

    jetAnalyzer->GetTracks(pisaRecoNode, all_track_list);

    std::vector<particlePair> track_pair_list_0;
    track_pair_list_0.clear();
    std::vector<particlePair> track_pair_list_1;
    track_pair_list_1.clear();
    std::vector<particlePair> track_pair_list_2;
    track_pair_list_2.clear();
    std::vector<particlePair> track_pair_list_3;
    track_pair_list_3.clear();

    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_constituent_list.size(); t++)
        {
            idTrue++;
            if(true_jet_constituent_list[t].charge == 0)
                {
                    continue;
                }
            float tPt        = true_jet_constituent_list[t].pT;
            float tEt        = true_jet_constituent_list[t].eT;
            float tEta       = true_jet_constituent_list[t].eta;
            float tPhi       = true_jet_constituent_list[t].phi;

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < all_track_list.size(); g++)
                {
                    idGood++;
                    if((true_jet_constituent_list[t].arm != all_track_list[g].arm) ||
		       (true_jet_constituent_list[t].charge != all_track_list[g].charge))
                        {
                            continue;
                        }
                    float gPt        = all_track_list[g].pT;
                    float gEt        = all_track_list[g].eT;
                    float gEta       = all_track_list[g].eta;
                    float gPhi       = all_track_list[g].phi;

                    float deltaR = jetAnalyzer->dR(tEta, gEta, tPhi, gPhi);

                    particlePair tempPair;
                    tempPair.id.first  = idTrue;
                    tempPair.id.second = idGood;
                    tempPair.pT.first  = tPt;
                    tempPair.pT.second = gPt;
                    tempPair.eT.first  = tEt;
                    tempPair.eT.second = gEt;
                    tempPair.deltaR    = deltaR;

                    bool passQuality     = all_track_list[g].passQuality;
                    bool passMatching    = all_track_list[g].passMatching;

                    float emcsdphi       = all_track_list[t].emcsdphi;
                    float emcsdz         = all_track_list[t].emcsdz;
                    float mom            = all_track_list[t].mom;
                    float energy         = all_track_list[t].energy;
                    int n0               = all_track_list[t].n0;

                    bool conversionEdge = jetAnalyzer->inEdge(gPhi);

                    bool conversionElectron = gPt < 4.5 && n0 >= 2 && (energy / mom) < 0.6;

                    bool emcMatching = sqrt((emcsdphi * emcsdphi) + (emcsdz * emcsdz)) < 3.0;
                    bool conversionEcore = emcMatching && energy < 0.2;

                    bool passConversions = !conversionEdge && !conversionElectron && !conversionEcore;

                    bool pass1 = passQuality;
                    bool pass2 = passQuality && passMatching;
                    bool pass3 = passQuality && passMatching && passConversions;


                    track_pair_list_0.push_back(tempPair);

                    if(pass1)
                        {
                            track_pair_list_1.push_back(tempPair);
                        }
                    if(pass2)
                        {
                            track_pair_list_2.push_back(tempPair);
                        }
                    if(pass3)
                        {
                            track_pair_list_3.push_back(tempPair);
                        }
                }
        }

    //Sort pair by ascending order of delta and do one-to-one matchin
    std::sort(track_pair_list_0.begin(), track_pair_list_0.end(), sortPair());

    std::vector<particlePair> unique_track_pair_list_0;
    unique_track_pair_list_0.clear();
    while (track_pair_list_0.size())
        {
            unique_track_pair_list_0.push_back(track_pair_list_0.front());
            track_pair_list_0.erase(std::remove_if(track_pair_list_0.begin(), track_pair_list_0.end(),
                                                   removePairId(track_pair_list_0.front())), track_pair_list_0.end());
        }

    for(unsigned int u = 0; u < unique_track_pair_list_0.size(); u++)
        {
            float truePt = unique_track_pair_list_0[u].pT.first;

            float deltaR = unique_track_pair_list_0[u].deltaR;
            hTDistance[0]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hTEfficiency[0]->Fill(truePt);
                }
        }

    /////////////////////////////////////////////////////////////////////////
    std::sort(track_pair_list_1.begin(), track_pair_list_1.end(), sortPair());

    std::vector<particlePair> unique_track_pair_list_1;
    unique_track_pair_list_1.clear();
    while (track_pair_list_1.size())
        {
            unique_track_pair_list_1.push_back(track_pair_list_1.front());
            track_pair_list_1.erase(std::remove_if(track_pair_list_1.begin(), track_pair_list_1.end(),
                                                   removePairId(track_pair_list_1.front())), track_pair_list_1.end());
        }

    for(unsigned int u = 0; u < unique_track_pair_list_1.size(); u++)
        {
            float truePt = unique_track_pair_list_1[u].pT.first;

            float deltaR = unique_track_pair_list_1[u].deltaR;
            hTDistance[1]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hTEfficiency[1]->Fill(truePt);
                }
        }

    /////////////////////////////////////////////////////////////////////////
    std::sort(track_pair_list_2.begin(), track_pair_list_2.end(), sortPair());

    std::vector<particlePair> unique_track_pair_list_2;
    unique_track_pair_list_2.clear();
    while (track_pair_list_2.size())
        {
            unique_track_pair_list_2.push_back(track_pair_list_2.front());
            track_pair_list_2.erase(std::remove_if(track_pair_list_2.begin(), track_pair_list_2.end(),
                                                   removePairId(track_pair_list_2.front())), track_pair_list_2.end());
        }

    for(unsigned int u = 0; u < unique_track_pair_list_2.size(); u++)
        {
            float truePt = unique_track_pair_list_2[u].pT.first;

            float deltaR = unique_track_pair_list_2[u].deltaR;
            hTDistance[2]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hTEfficiency[2]->Fill(truePt);
                }
        }

    /////////////////////////////////////////////////////////////////////////
    std::sort(track_pair_list_3.begin(), track_pair_list_3.end(), sortPair());

    std::vector<particlePair> unique_track_pair_list_3;
    unique_track_pair_list_3.clear();
    while (track_pair_list_3.size())
        {
            unique_track_pair_list_3.push_back(track_pair_list_3.front());
            track_pair_list_3.erase(std::remove_if(track_pair_list_3.begin(), track_pair_list_3.end(),
                                                   removePairId(track_pair_list_3.front())), track_pair_list_3.end());
        }

    for(unsigned int u = 0; u < unique_track_pair_list_3.size(); u++)
        {
            float truePt = unique_track_pair_list_3[u].pT.first;

            float deltaR = unique_track_pair_list_3[u].deltaR;
            hTDistance[3]->Fill(deltaR);

            if(deltaR < minChargedDeltaR)
                {
                    hTEfficiency[3]->Fill(truePt);
                }
        }
}


void JetConstituentPerformance::fillHistograms()
{
    for(int t = 0; t < pythia_true_jet_list.size(); t++)
        {
            float pT = pythia_true_jet_list[t].pT;
            hPtTrueJet->Fill(pT);
        }

    for(int c = 0; c < pythia_true_jet_constituent_list.size(); c++)
        {
            float jetPt = pythia_true_jet_constituent_list[c].jetPt;
            float pT = pythia_true_jet_constituent_list[c].pT;
            int kf = pythia_true_jet_constituent_list[c].id;

            if((kf == PY_PI) || (kf == -PY_PI))
                {
                    hPtPion->Fill(jetPt, pT / jetPt);
                }
            if(kf == PY_GAMMA)
                {
                    hPtPhoton->Fill(jetPt, pT / jetPt);
                }
            if((kf == PY_K) || (kf == -PY_K))
                {
                    hPtKaon->Fill(jetPt, pT / jetPt);
                }
            if((kf == PY_P) || (kf == -PY_P))
                {
                    hPtProton->Fill(jetPt, pT / jetPt);
                }
            if((kf == PY_ELECTRON) || (kf == -PY_ELECTRON))
                {
                    hPtElectron->Fill(jetPt, pT / jetPt);
                }
            if((kf == 2112) || (kf == -2112))
                {
                    hPtNeutron->Fill(jetPt, pT / jetPt);
                }
            if((kf == 130) || (kf == -130))
                {
                    hPtK0L->Fill(jetPt, pT / jetPt);
                }
        }
}

void JetConstituentPerformance::fillTrees(bool what)
{
    if(what)
        {
            tPythiaParticle->Fill();
            tPythiaTrueJet->Fill();
            tPythiaTrueJetConstituent->Fill();

            tPisaRecoAllTrack->Fill();
            tPisaRecoAllCluster->Fill();

            tPisaRecoGoodTrack->Fill();
            tPisaRecoGoodCluster->Fill();
            tPisaRecoGoodParticle->Fill();
        }
}

int JetConstituentPerformance::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetConstituentPerformance::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nTotalEvents);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nTotalEvents << endl;

    jetAnalyzer->MyStatistics();

    if(nTotalEvents > 0)
        {
            cout << "This is the end. If everthing was ok, you will see this message: " << endl;
            cout << "123ALLDONE321" << endl;
        }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}















































