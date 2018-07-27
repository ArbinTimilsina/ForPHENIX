//General PHENIX tools
#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <getClass.h>

//HepMc tools
#include <PHHepMCGenEvent.h>

//For vertex during sim reco
#include <VtxOut.h>

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

//My source file
#include "FakeJetStudy.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
FakeJetStudy::FakeJetStudy()
    : SubsysReco("FakeJetStudy")
{
    nTotalEvents = 0;
    nTotalEventsCentrality1 = 0;
    nTotalEventsCentrality2 = 0;
    nTotalEventsCentrality3 = 0;
    nTotalEventsCentrality4 = 0;

    nTrueJetEvents = 0;
    nTrueJetEventsCentrality1 = 0;
    nTrueJetEventsCentrality2 = 0;
    nTrueJetEventsCentrality3 = 0;
    nTrueJetEventsCentrality4 = 0;

    nNoJetEvents = 0;
    nNoJetEventsCentrality1 = 0;
    nNoJetEventsCentrality2 = 0;
    nNoJetEventsCentrality3 = 0;
    nNoJetEventsCentrality4 = 0;

    treesFill = false;

    return;
}

int FakeJetStudy::Init(PHCompositeNode *topNode)
{
    //To access various nodes
    se = Fun4AllServer::instance();
    rc = recoConsts::instance();

    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(false);
    jetAnalyzer->SetCuAu(true);
    jetAnalyzer->MyInit();

    //Output file name
    outfile = new TFile("FakeJetStudy.root", "RECREATE");

    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    hEvents = new TH1F("hEvents", "Number of events", 15, 0, 15);
    hImpactParameter = new TH1F("hImpactParameter", "sHIJING Impact Parameter", 210, -1, 20);
    hEventPlane = new TH1F("hEventPlane", "sHIJING Event Plane", 140, -0.5, 6.5);

    for (unsigned int c = 0; c <= 4; c++)
        {
            hTracksZed[c] = new TH1F(Form("hTracksZed_%u", c),
                                     "#frac{dN}{dZed} of tracks going into Jet reconstruction", 400, -100.0, 100.0);
            hTracksPhi[c] = new TH1F(Form("hTracksPhi_%u", c),
                                     "#frac{dN}{d#phi} of tracks going into Jet reconstruction", 200, -1.0, 4.0);

            hClustersEta[c] = new TH1F(Form("hClustersEta_%u", c),
                                       "#frac{dN}{d#eta} of clusters going into Jet reconstruction", 160, -0.4, 0.4);
            hClustersPhi[c] = new TH1F(Form("hClustersPhi_%u", c),
                                       "#frac{dN}{d#phi} of clusters going into Jet reconstruction", 200, -1.0, 4.0);
            hClustersTower[c] = new TH1F(Form("hClustersTower_%u", c),
                                         "Hits vs Tower-id", NTOWER, 0, NTOWER);

            hDistance[c] = new TH1F(Form("hDistance_%u", c), "#DeltaR between True Jets and Reco/Embedded Jets",
                                    400, 0.0, 0.4);
            hDistanceVsPtTrue[c] = new TH2F(Form("hDistanceVsPtTrue_%u", c), "#DeltaR vs p_{T, true} for Reco/Embedded Jets",
                                            40, 10, 50, 40, 0.0, 0.4);

            hPtTrue[c] = new TH1F(Form("hPtTrue_%u", c),
                                  "p_{T} of True Jets", NPTBINS_TRUE, PTBINS_TRUE);
            hPtTrueMatched[c] = new TH1F(Form("hPtTrueMatched_%u", c),
                                         "p_{T} of Matched True Jets", NPTBINS_TRUE, PTBINS_TRUE);
            hPtReco[c] = new TH1F(Form("hPtReco_%u", c),
                                  "p_{T} of Reco Jets", NPTBINS_RECO, PTBINS_RECO);
            hPtRecoMatched[c] = new TH1F(Form("hPtRecoMatched_%u", c),
                                         "p_{T} of Matched Reco Jets", NPTBINS_RECO, PTBINS_RECO);
            hPtRecoNotMatched[c] = new TH1F(Form("hPtRecoNotMatched_%u", c),
                                            "p_{T} of Not Matched Reco Jets", NPTBINS_RECO, PTBINS_RECO);

            hResponseMatrix[c] = new TH2F(Form("hResponseMatrix_%u", c), "p_{T, true} vs p_{T, reco} for Matched Jets",
                                          NPTBINS_RECO, PTBINS_RECO, NPTBINS_TRUE, PTBINS_TRUE);

            hPtFake[c] = new TH1F(Form("hPtFake_%u", c),
                                  "p_{T} of \"Fake\" Jets", NPTBINS_RECO, PTBINS_RECO);
        }

    //*****************************************************************************************************************************
    //Trees
    //*****************************************************************************************************************************
    if(treesFill)
        {
            tTrueJets = new TTree("tTrueJets", "True Jets from sHIJING");
            tTrueJets->Branch("tTrueJets", &true_jets);

            tRecoAllTracks = new TTree("tRecoAllTracks", "Simulated All Tracks");
            tRecoAllTracks->Branch("tRecoAllTracks", &reco_all_tracks);

            tRecoAllClusters = new TTree("tRecoAllClusters", "Simulated All Clusters");
            tRecoAllClusters->Branch("tRecoAllClusters", &reco_all_clusters);

            tRecoGoodTracks = new TTree("tRecoGoodTracks", "Simulated Tracks with all cuts");
            tRecoGoodTracks->Branch("tRecoGoodTracks", &reco_good_tracks);

            tRecoGoodClusters = new TTree("tRecoGoodClusters", "Simulated Clusters with all cuts");
            tRecoGoodClusters->Branch("tRecoGoodClusters", &reco_good_clusters);

            tRecoGoodParticles = new TTree("tRecoGoodParticles", "Simulated Tracks and Clusters with all cuts");
            tRecoGoodParticles->Branch("tRecoGoodParticles", &reco_good_particles);

            tRecoJets = new TTree("tRecoJets", "Reconstrctred Jets from simulated Tracks and Clusters with all cuts");
            tRecoJets->Branch("tRecoJets", &reco_jets);

            tMatchedJets = new TTree("tMatchedJets", "Matched Reco Jets");
            tMatchedJets->Branch("tMatchedJets", &matched_jets);

            tNotMatchedJets = new TTree("tNotMatchedJets", "Not Matched Reco Jets");
            tNotMatchedJets->Branch("tNotMatchedJets", &not_matched_jets);

            //Fake jet
            tShuffledParticles = new TTree("tShuffledParticles", "Shuffled particles for fake jet");
            tShuffledParticles->Branch("tShuffledParticles", &shuffled_particles);

            tFakeJets = new TTree("tFakeJets", "Fake Jet yield");
            tFakeJets->Branch("tFakeJets", &fake_jets);

            tFakeConstituents = new TTree("tFakeConstituents", "Fake Jet constituent");
            tFakeConstituents->Branch("tFakeConstituents", &fake_jet_constituents);
        }

    return 0;
}

int FakeJetStudy::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int FakeJetStudy::ResetEvent(PHCompositeNode *topNode)
{
    true_jets.clear();

    reco_all_tracks.clear();
    reco_all_clusters.clear();

    reco_good_tracks.clear();
    reco_good_clusters.clear();
    reco_good_particles.clear();
    reco_jets.clear();
    reco_jets_for_fake.clear();

    matched_jets.clear();
    not_matched_jets.clear();

    shuffled_particles.clear();
    fake_jets.clear();
    fake_jet_constituents.clear();

    return EVENT_OK;
}


int FakeJetStudy::process_event(PHCompositeNode *topNode)
{
    //Find the HepMC data
    PHHepMCGenEvent *genEvt = findNode::getClass<PHHepMCGenEvent>(topNode, "PHHepMCGenEvent");
    HepMC::GenEvent *sHijing =  genEvt->getEvent();
    if (!sHijing)
        {
            cout << "No sHijing! No sense continuing" << endl;
            exit(1);
        }

    //Global things
    HepMC::HeavyIon* heavyIon = sHijing->heavy_ion();

    //Impact parameter in fm
    float impactParameter = (float)heavyIon->impact_parameter();
    hImpactParameter->Fill(impactParameter);

    //Calculate centrality
    centralityBin = getCentralityBin(impactParameter);

    //Azimuthal angle of event plane
    float eventPlane = (float)heavyIon->event_plane_angle();
    hEventPlane->Fill(eventPlane);

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

    getTrueJet(sHijing,
               true_jets, hPtTrue);

    if(true_jets.size() != 0)
        {
            nTrueJetEvents++;
            if(centralityBin == 1)
                {
                    nTrueJetEventsCentrality1++;
                }
            if(centralityBin == 2)
                {
                    nTrueJetEventsCentrality2++;
                }
            if(centralityBin == 3)
                {
                    nTrueJetEventsCentrality3++;
                }
            if(centralityBin == 4)
                {
                    nTrueJetEventsCentrality4++;
                }
        }

    getRecoJet(topNode,
               reco_all_tracks, reco_all_clusters,
               reco_good_tracks, reco_good_clusters, reco_good_particles,
               reco_jets, reco_jets_for_fake, hPtReco);

    getMatchedJet(true_jets, reco_jets,
                  matched_jets, not_matched_jets,
                  hDistance, hDistanceVsPtTrue, hPtTrueMatched, hPtRecoMatched, hPtRecoNotMatched, hResponseMatrix);

    //Get Fake jets in non-jet events (pT>10.45), Anti-kt, R=0.2
    if(reco_jets_for_fake.size() == 0)
        {
            jetAnalyzer->GetFakeJets(reco_good_particles, 0.2, 3.0, minPt, 0.2, 0.7,
                                     shuffled_particles, fake_jets, fake_jet_constituents);

            nNoJetEvents++;
            if(centralityBin == 1)
                {
                    nNoJetEventsCentrality1++;
                }
            if(centralityBin == 2)
                {
                    nNoJetEventsCentrality2++;
                }
            if(centralityBin == 3)
                {
                    nNoJetEventsCentrality3++;
                }
            if(centralityBin == 4)
                {
                    nNoJetEventsCentrality4++;
                }
        }

    fillTrees(treesFill);

    ///////////////////////////////////////////////////////////////////////////////////////
    //Fill histograms
    ///////////////////////////////////////////////////////////////////////////////////////
    for (unsigned int t = 0; t < reco_good_tracks.size(); t++)
        {
            float zedDC           = reco_good_tracks[t].zedDC;
            float phiDC           = reco_good_tracks[t].phiDC;

            hTracksZed[0]->Fill(zedDC);
            hTracksPhi[0]->Fill(phiDC);

            if(centralityBin == 1)
                {
                    hTracksZed[1]->Fill(zedDC);
                    hTracksPhi[1]->Fill(phiDC);
                }
            if(centralityBin == 2)
                {
                    hTracksZed[2]->Fill(zedDC);
                    hTracksPhi[2]->Fill(phiDC);
                }
            if(centralityBin == 3)
                {
                    hTracksZed[3]->Fill(zedDC);
                    hTracksPhi[3]->Fill(phiDC);
                }
            if(centralityBin == 4)
                {
                    hTracksZed[4]->Fill(zedDC);
                    hTracksPhi[4]->Fill(phiDC);
                }
        }

    for (unsigned int c = 0; c < reco_good_clusters.size(); c++)
        {
            float eta = reco_good_clusters[c].eta;
            float phi = reco_good_clusters[c].phi;
            int towerId = reco_good_clusters[c].towerId;

            hClustersEta[0]->Fill(eta);
            hClustersPhi[0]->Fill(phi);
            hClustersTower[0]->Fill(towerId);

            if(centralityBin == 1)
                {
                    hClustersEta[1]->Fill(eta);
                    hClustersPhi[1]->Fill(phi);
                    hClustersTower[1]->Fill(towerId);
                }
            if(centralityBin == 2)
                {
                    hClustersEta[2]->Fill(eta);
                    hClustersPhi[2]->Fill(phi);
                    hClustersTower[2]->Fill(towerId);
                }
            if(centralityBin == 3)
                {
                    hClustersEta[3]->Fill(eta);
                    hClustersPhi[3]->Fill(phi);
                    hClustersTower[3]->Fill(towerId);
                }
            if(centralityBin == 4)
                {
                    hClustersEta[4]->Fill(eta);
                    hClustersPhi[4]->Fill(phi);
                    hClustersTower[4]->Fill(towerId);
                }
        }

    for (unsigned int f = 0; f < fake_jets.size(); f++)
        {
            float pT  = fake_jets[f].pT;

            hPtFake[0]->Fill(pT);

            if(centralityBin == 1)
                {
                    hPtFake[1]->Fill(pT);
                }
            if(centralityBin == 2)
                {
                    hPtFake[2]->Fill(pT);
                }
            if(centralityBin == 3)
                {
                    hPtFake[3]->Fill(pT);
                }
            if(centralityBin == 4)
                {
                    hPtFake[4]->Fill(pT);
                }
        }
    return 0;
}

void FakeJetStudy::getTrueJet(HepMC::GenEvent *sHijingEvt,
                              std::vector<jets>& true_jets_list,
                              TH1F *hTrue[5])
{
    if(nTotalEvents == 5)
        {
            cout << endl;
            cout << "***********************************************************************" << endl;
            cout << "sHIJING particle information:" << endl;
            cout << "***********************************************************************" << endl;
            cout << "Particle's:" << setw(15) << "pT" << setw(15) << "eta" << setw(15) << "phi" << endl;
        }

    for ( HepMC::GenEvent::particle_iterator p = sHijingEvt->particles_begin(); p != sHijingEvt->particles_end(); ++p )
        {
            const HepMC::FourVector& momVector = (*p)->momentum();

            float rho = sqrt(pow(momVector.px(), 2) + pow(momVector.py(), 2) + pow(momVector.pz(), 2));
            float jetEta = 0.5 * log((rho + momVector.pz()) / (rho - momVector.pz()));

            float jetPt = sqrt(pow(momVector.px(), 2) + pow(momVector.py(), 2));
            float jetPhi = jetAnalyzer->phiReduce(atan2(momVector.py(), momVector.px()));

            int arm = 1;
            if (jetPhi > 1.57)
                {
                    arm = 0;
                }

            jets tempJet;
            tempJet.arm           = arm;
            tempJet.centralityBin = centralityBin;
            tempJet.pT            = jetPt;
            tempJet.eta           = jetEta;
            tempJet.phi           = jetPhi;
            tempJet.nc            = -99.9;
            tempJet.cf            = -99.9;
            tempJet.nf            = -99.9;
            tempJet.disc          = -99.9;

            int status = (*p)->status();
            int id = (*p)->pdg_id();

            if((status == 103) && (id == 2000000))
                {
                    bool acceptanceCut = meetPhenixAcceptance(jetEta, jetPhi);
                    if (acceptanceCut && (jetPt > MINPT_TRUE))
                        {
                            true_jets_list.push_back(tempJet);

                            hTrue[0]->Fill(tempJet.pT);

                            if(centralityBin == 1)
                                {
                                    hTrue[1]->Fill(tempJet.pT);
                                }
                            if(centralityBin == 2)
                                {
                                    hTrue[2]->Fill(tempJet.pT);
                                }
                            if(centralityBin == 3)
                                {
                                    hTrue[3]->Fill(tempJet.pT);
                                }
                            if(centralityBin == 4)
                                {
                                    hTrue[4]->Fill(tempJet.pT);
                                }
                        }
                }
        }
}

void FakeJetStudy::getRecoJet(PHCompositeNode *topNode,
                              std::vector<tracks>& all_track_list,
                              std::vector<clusters>& all_cluster_list,
                              std::vector<tracks>& good_track_list,
                              std::vector<clusters>& good_cluster_list,
                              std::vector<particles>& good_particle_list,
                              std::vector<jets>& reco_jet_list,
                              std::vector<jets>& reco_jet_for_fake_list,
                              TH1F *hReco[5])
{
    PHCompositeNode *pisaRecoNode = se->topNode(rc->get_CharFlag("PISARECO_TOPNODE"));

    //For ZVertex used in reconstruction
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
    jetAnalyzer->GetClusters(pisaRecoNode, all_cluster_list, zvertex);

    jetAnalyzer->GetParticles(all_track_list, all_cluster_list, good_track_list, good_cluster_list, good_particle_list);

    std::vector<particles> constituent_list;
    constituent_list.clear();

    jetAnalyzer->GetAntiKtCommon(good_particle_list, R, nc, minPt, minCf, maxCf, reco_jet_list, constituent_list);

    for(unsigned int r = 0; r < reco_jet_list.size(); r++)
        {
            float pT = reco_jet_list[r].pT;

            if(pT > MINPT_RECO)
                {
                    reco_jet_for_fake_list.push_back(reco_jet_list[r]);
                }

            hReco[0]->Fill(pT);

            if(centralityBin == 1)
                {
                    hReco[1]->Fill(pT);
                }
            if(centralityBin == 2)
                {
                    hReco[2]->Fill(pT);
                }
            if(centralityBin == 3)
                {
                    hReco[3]->Fill(pT);
                }
            if(centralityBin == 4)
                {
                    hReco[4]->Fill(pT);
                }
        }
}

void FakeJetStudy::getMatchedJet(std::vector<jets> true_jet_list,
                                 std::vector<jets> reco_jet_list,
                                 std::vector<jets>& matched_jet_list,
                                 std::vector<jets>& not_matched_jet_list,
                                 TH1F *hDistance[5],
                                 TH2F *hDistanceVsPtTrue[5],
                                 TH1F *hTrueMatched[5],
                                 TH1F *hRecoMatched[5],
                                 TH1F *hRecoNotMatched[5],
                                 TH2F *hMatrix[5])
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

            unsigned int idReco = 0;
            for(unsigned int r = 0; r < reco_jet_list.size(); r++)
                {
                    idReco++;
                    float rArm       = reco_jet_list[r].arm;
                    float rPt        = reco_jet_list[r].pT;
                    float rEta       = reco_jet_list[r].eta;
                    float rPhi       = reco_jet_list[r].phi;
                    float rNc        = reco_jet_list[r].nc;
                    float rCf        = reco_jet_list[r].cf;
                    float rNf        = reco_jet_list[r].nf;
                    float rDisc      = reco_jet_list[r].disc;

                    if (tArm == rArm)
                        {
                            float deltaR = jetAnalyzer->dR(tEta, rEta, tPhi, rPhi);

                            jetPair tempPair;
                            tempPair.id.first    = idTrue;
                            tempPair.id.second   = idReco;
                            tempPair.pT.first    = tPt;
                            tempPair.pT.second   = rPt;
                            tempPair.eta.first   = tEta;
                            tempPair.eta.second  = rEta;
                            tempPair.phi.first   = tPhi;
                            tempPair.phi.second  = rPhi;
                            tempPair.nc.first    = tNc;
                            tempPair.nc.second   = rNc;
                            tempPair.cf.first    = tCf;
                            tempPair.cf.second   = rCf;
                            tempPair.nf.first    = tNf;
                            tempPair.nf.second   = rNf;
                            tempPair.disc.first  = tDisc;
                            tempPair.disc.second = rDisc;
                            tempPair.arm         = tArm;
                            tempPair.deltaR      = deltaR;

                            jet_pair_list.push_back(tempPair);
                        }
                }
        }

    //Sort pair by ascending order of deltaR
    std::sort(jet_pair_list.begin(), jet_pair_list.end(), sortPair());

    //Matching: 1 true jet can be matched to 2 reco jets (due to splitting)- and save as unique pair
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

            if (deltaR < minDeltaR)
                {
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

                    hTrueMatched[0]->Fill(truePt);
                    hRecoMatched[0]->Fill(recoPt);
                    hMatrix[0]->Fill(recoPt, truePt);

                    if(centralityBin == 1)
                        {
                            hTrueMatched[1]->Fill(truePt);
                            hRecoMatched[1]->Fill(recoPt);
                            hMatrix[1]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 2)
                        {
                            hTrueMatched[2]->Fill(truePt);
                            hRecoMatched[2]->Fill(recoPt);
                            hMatrix[2]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 3)
                        {
                            hTrueMatched[3]->Fill(truePt);
                            hRecoMatched[3]->Fill(recoPt);
                            hMatrix[3]->Fill(recoPt, truePt);
                        }
                    if(centralityBin == 4)
                        {
                            hTrueMatched[4]->Fill(truePt);
                            hRecoMatched[4]->Fill(recoPt);
                            hMatrix[4]->Fill(recoPt, truePt);
                        }
                }
        }

    for(unsigned int r = 0; r < reco_jet_list.size(); r++)
        {
            jets tempJet;
            tempJet.arm           = reco_jet_list[r].arm;
            tempJet.centralityBin = centralityBin;
            tempJet.pT            = reco_jet_list[r].pT;
            tempJet.eta           = reco_jet_list[r].eta;
            tempJet.phi           = reco_jet_list[r].phi;
            tempJet.nc            = reco_jet_list[r].nc;
            tempJet.cf            = reco_jet_list[r].cf;
            tempJet.nf            = reco_jet_list[r].nf;
            tempJet.disc          = reco_jet_list[r].disc;

            bool match = false;
            for(unsigned int m = 0; m < matched_jet_list.size(); m++)
                {
                    float mArm       = matched_jet_list[m].arm;
                    float mPt        = matched_jet_list[m].pT;
                    float mEta       = matched_jet_list[m].eta;
                    float mPhi       = matched_jet_list[m].phi;
                    float mNc        = matched_jet_list[m].nc;
                    float mCf        = matched_jet_list[m].cf;
                    float mNf        = matched_jet_list[m].nf;
                    float mDisc      = matched_jet_list[m].disc;

                    bool check = (mArm == tempJet.arm) &&
			(fabs(mPt - tempJet.pT) < 0.0001) && (fabs(mEta - tempJet.eta) < 0.0001) &&
			(fabs(mPhi - tempJet.phi) < 0.0001) && (fabs(mNc - tempJet.nc) < 0.0001) &&
			(fabs(mCf - tempJet.cf) < 0.0001) && (fabs(mNf - tempJet.nf) < 0.0001) &&
			(fabs(mDisc - tempJet.disc) < 0.0001);

                    if (check)
                        {
                            match = true;
                        }
                }

            if(!match)
                {
                    not_matched_jet_list.push_back(tempJet);

                    hRecoNotMatched[0]->Fill(tempJet.pT);

                    if(centralityBin == 1)
                        {
                            hRecoNotMatched[1]->Fill(tempJet.pT);
                        }
                    if(centralityBin == 2)
                        {
                            hRecoNotMatched[2]->Fill(tempJet.pT);
                        }
                    if(centralityBin == 3)
                        {
                            hRecoNotMatched[3]->Fill(tempJet.pT);
                        }
                    if(centralityBin == 4)
                        {
                            hRecoNotMatched[4]->Fill(tempJet.pT);
                        }
                }
        }
}


void FakeJetStudy::fillTrees(bool what)
{
    if(what)
        {
            tTrueJets->Fill();

            tRecoAllTracks->Fill();
            tRecoAllClusters->Fill();

            tRecoGoodTracks->Fill();
            tRecoGoodClusters->Fill();
            tRecoGoodParticles->Fill();
            tRecoJets->Fill();

            tMatchedJets->Fill();
            tNotMatchedJets->Fill();

            tShuffledParticles->Fill();
            tFakeJets->Fill();
            tFakeConstituents->Fill();
        }
}

int FakeJetStudy::End(PHCompositeNode *topNode)
{
    //Bins 1 and 2
    hEvents->SetBinContent(1, nTotalEvents);
    hEvents->SetBinContent(2, nTrueJetEvents);

    //Bins 4 and 5
    hEvents->SetBinContent(4, nTotalEventsCentrality1);
    hEvents->SetBinContent(5, nTrueJetEventsCentrality1);

    //Bins 7 and 8
    hEvents->SetBinContent(7, nTotalEventsCentrality2);
    hEvents->SetBinContent(8, nTrueJetEventsCentrality2);

    //Bins 10 and 11
    hEvents->SetBinContent(10, nTotalEventsCentrality3);
    hEvents->SetBinContent(11, nTrueJetEventsCentrality3);

    //Bins 13 and 14
    hEvents->SetBinContent(13, nTotalEventsCentrality4);
    hEvents->SetBinContent(14, nTrueJetEventsCentrality4);

    //Bin 3, 6, 9, 12 and 15 for Fakes
    hEvents->SetBinContent(3, nNoJetEvents);
    hEvents->SetBinContent(6, nNoJetEventsCentrality1);
    hEvents->SetBinContent(9, nNoJetEventsCentrality2);
    hEvents->SetBinContent(12, nNoJetEventsCentrality3);
    hEvents->SetBinContent(15, nNoJetEventsCentrality4);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nTotalEvents << endl;
    cout << "Total events with true jet = " << nTrueJetEvents << endl;

    if(nTotalEvents > 0)
        {
            cout << "This is the end. If everthing was ok, you will see this message: " << endl;
            cout << "123ALLDONE321" << endl;
        }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    outfile->Write();
    outfile->Close();

    return 0;
}
