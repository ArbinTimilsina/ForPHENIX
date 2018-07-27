//General PHENIX tools
#include <getClass.h>
#include <PHCompositeNode.h>
#include <phool.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//My source file
#include "JetAnalyzer.h"

//For Quality mask
#include "TrackQualityCuAu.h"
#include "TrackQualityPP.h"

//Used for PC3 matching
#include "Pc3Matching.h"

//Used for EMCal matching
#include "EmcMatching.h"

//For bitmask
#define X1_USED 1
#define X2_USED 2
#define UV_FOUND_UNIQUE 12

//For MB trigger
static const int MB_PP_NARROWVTX   = 0x00000010;
static const int MB_CUAU_NARROWVTX = 0x00000004;

//For ERT trigger
static const int ERT_BIT_MASK_C = 0x00000100;

using namespace std;
using namespace findNode;
using namespace fastjet;
//================================ Constructor ================================
//Here we can initiate some variables
JetAnalyzer::JetAnalyzer(std::string outfilename)
    : SubsysReco("JetAnalyzer"),
      verbo(1),
      outfname(outfilename),
      writeTree(false)
{
    outfile = new TFile(outfname.c_str(), "RECREATE");

    for (unsigned int c = 0; c < 5; c++)
        {
            nTotalEvents[c] = 0;
            nGoodEvents[c] = 0;
            nNoJetEvents[c] = 0;
            nNoJetEventsEast[c] = 0;
            nNoJetEventsWest[c] = 0;
            nNoJetEventsFidTight[c] = 0;
            nNoJetEventsNc[c] = 0;
            nNoJetEventsCf[c] = 0;
            nNoJetEventsNcCf[c] = 0;
            nNoJetEventsTrClTight[c] = 0;

            nJets[c] = 0;
            nJetsHighPt[c] = 0;

            nVertexEvents[c] = 0;
            nTotalZvertex[c] = 0.0;

            nNoR3JetEvents[c] = 0;
        }

    nEffectiveEvents = 0.0;

    nAllTracks = 0;
    nAllClusters = 0;
    nChargedParticles = 0;
    nNeutralParticles = 0;

    nPassPt = 0;
    nPassQuality = 0;
    nPassEMCMatching = 0;
    nPassPC3Matching = 0;
    nPassMatching = 0;

    nPassFirst = 0;
    nPassDcPcGhost = 0;
    nPassDcPcConversion = 0;
    nPassPairCut = 0;
    nPassConversionEdge = 0;
    nPassConversionElectron = 0;
    nPassConversionEcore = 0;
    nPassConversion = 0;

    nPassIsValid = 0;
    nPassNotHot = 0;
    nPassNotDead = 0;
    nPassNotUncalib = 0;
    nPassNotBad = 0;
    nPassEnergy = 0;
    nPassEverything = 0;

    nCentralityEvents = 0;
    nTotalCentrality = 0.0;
    nTotalCentralityAfterVertex = 0.0;

    //Jet algorithm definition
    //Anit-kt
    antikt_00 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.15, fastjet::E_scheme, fastjet::Best);
    antikt_01 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.2, fastjet::E_scheme, fastjet::Best);
    antikt_02 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.25, fastjet::E_scheme, fastjet::Best);
    antikt_03 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.3, fastjet::E_scheme, fastjet::Best);

    InitTrees(writeTree);
    InitHistograms();

    return;
}

int JetAnalyzer::Init(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    //Read warnmap, quality and matching stuff
    MyInit();

    return 0;
}

int JetAnalyzer::InitRun(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>> InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    RunHeader *run_header = getClass<RunHeader> (topNode, "RunHeader");
    if (!run_header)
        {
            cout << "No RunHeader! No sense continuing" << endl;
            exit(1);
        }

    runNumber = run_header->get_RunNumber();

    triggerHelper = new TriggerHelper(topNode);

    return EVENT_OK;
}


int JetAnalyzer::process_event(PHCompositeNode *topNode)
{
    //Bad run list- this needs to happen here (in process event)
    if(isCuAu)
        {
            //CuAu
            //First run was 372402 and last run was 377310
            if(runNumber < 372402 || runNumber > 377310)
                {
                    return DISCARDEVENT;
                }

            //CuAu bad run list: 30 Runs out of 451
            if(runNumber == 372402 || runNumber == 372524 || runNumber == 372525 || runNumber == 372531 || runNumber == 372533 ||
	       runNumber == 372536 || runNumber == 372647 || runNumber == 372648 || runNumber == 372959 || runNumber == 372961 ||

	       runNumber == 373407 || runNumber == 373655 || runNumber == 373672 ||

	       runNumber == 374428 ||

	       runNumber == 375773 || runNumber == 375774 || runNumber == 375906 || runNumber == 375953 || runNumber == 375957 ||

	       runNumber == 376433 || runNumber == 376434 || runNumber == 376435 || runNumber == 376620 ||

	       runNumber == 377155 || runNumber == 377156 || runNumber == 377157 || runNumber == 377167 || runNumber == 377171 ||
	       runNumber == 377172 || runNumber == 377173)
                {
                    return DISCARDEVENT;
                }
        }
    else
        {
            //p+p 200 GeV
            //First run was 357665 and last run was 363228
            if(runNumber < 357665 || runNumber > 363228)
                {
                    return DISCARDEVENT;
                }

            //p+p bad run list: 82 runs out of 429 bad (pro99). For pro101, 50 +3 out of 328 bad
            if(runNumber == 358661 || runNumber == 358663 || runNumber == 358665 || runNumber == 358667 || runNumber == 358710 ||
	       runNumber == 358711 || runNumber == 358717 || runNumber == 358719 || runNumber == 358720 || runNumber == 358722 ||
	       runNumber == 358724 || runNumber == 358725 || runNumber == 358740 || runNumber == 358742 || runNumber == 358743 ||
	       runNumber == 358749 || runNumber == 358750 || runNumber == 358751 || runNumber == 358752 || runNumber == 358754 ||
	       runNumber == 358758 || runNumber == 358759 || runNumber == 358767 || runNumber == 358768 || runNumber == 358771 ||
	       runNumber == 358772 || runNumber == 358773 || runNumber == 358776 || runNumber == 358777 || runNumber == 358778 ||
	       runNumber == 358779 || runNumber == 358780 || runNumber == 358782 || runNumber == 358783 || runNumber == 358898 ||
	       runNumber == 358899 || runNumber == 358900 || runNumber == 358903 || runNumber == 358904 || runNumber == 358924 ||
	       runNumber == 358985 || runNumber == 358986 || runNumber == 358988 || runNumber == 358991 || runNumber == 358992 ||
	       runNumber == 358996 || runNumber == 358997 || runNumber == 358998 || runNumber == 359002 || runNumber == 359060 ||
	       runNumber == 359061 || runNumber == 359062 || runNumber == 359064 ||

	       runNumber == 359293 || runNumber == 359520 || runNumber == 359696 || runNumber == 359791 ||

	       runNumber == 360075 || runNumber == 360076 || runNumber == 360077 || runNumber == 360079 || runNumber == 360081 ||
	       runNumber == 360082 || runNumber == 360083 || runNumber == 360088 || runNumber == 360089 || runNumber == 360125 ||
	       runNumber == 360126 || runNumber == 360128 || runNumber == 360132 || runNumber == 360135 || runNumber == 360136 ||
	       runNumber == 360138 || runNumber == 360139 || runNumber == 360140 || runNumber == 360141 ||

	       runNumber == 360510 || runNumber == 361244 || runNumber == 361640 || runNumber == 361641 || runNumber == 362214 ||
	       runNumber == 362260 ||

	       //Post QM15- bad mean vertex
	       runNumber == 360475 || runNumber == 360501 || runNumber == 363196)
                {
                    return DISCARDEVENT;
                }
        }

    //Make trigger selection
    TrigLvl1 *lvl1trigger = getClass<TrigLvl1>(topNode, "TrigLvl1");

    //MB trigger
    if(isMB)
        {
            int mbFired = 0;
            if(isCuAu)
                {
                    if(lvl1trigger->get_lvl1_trigscaled()&MB_CUAU_NARROWVTX)
                        {
                            mbFired = 1;
                        }
                }
            else
                {
                    if(lvl1trigger->get_lvl1_trigscaled()&MB_PP_NARROWVTX)
                        {
                            mbFired = 1;
                        }

                }

            if(!mbFired)
                {
                    return DISCARDEVENT;
                }
        }
    else //ERT trigger for p+p
        {
            int ertFired = 0;
            if(lvl1trigger->get_lvl1_trigscaled()&ERT_BIT_MASK_C)
                {
                    ertFired = 1;
                }

            if(!ertFired)
                {
                    return DISCARDEVENT;
                }
        }

    //Get the data I need...
    phglobal                               = getClass<PHGlobal>                (topNode, "PHGlobal");
    if (nTotalEvents[0] == 1 && !phglobal)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    evtheader                              = getClass<EventHeader>             (topNode, "EventHeader");
    if (nTotalEvents[0] == 1 && !evtheader)
        {
            cout << "No EventHeader! No sense continuing" << endl;
            exit(1);
        }

    //Centrality
    centrality = phglobal->getCentrality();

    //Centrality cut- some events can have undefined centrality
    if(isCuAu && (centrality < 0.0 || centrality > 100.0))
        {
            return DISCARDEVENT;
        }

    hCentrality->Fill(centrality);
    nCentralityEvents++;
    nTotalCentrality += centrality;

    centralityBin = getCentralityBin(centrality);

    //Vertex
    zvertex = phglobal->getBbcZVertex();

    if(fabs(zvertex) < 100)
        {
            hVertex[0]->Fill(zvertex);
            nVertexEvents[0]++;
            nTotalZvertex[0] += zvertex;
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    hVertex[centralityBin]->Fill(zvertex);
                    nVertexEvents[centralityBin]++;
                    nTotalZvertex[centralityBin] += zvertex;
                }
        }

    //These are the events for the trigger selection and good runs
    nTotalEvents[0]++;
    if(centralityBin >= 1 && centralityBin <= 4)
        {
            nTotalEvents[centralityBin]++;
        }
    // Informational message...
    if (nTotalEvents[0] % 1000 == 0 && verbosity)
        {
            if (verbo)
                {
                    cout << "Events analyzed for run " << runNumber << " = " << nTotalEvents[0] << endl;
                }
        }

    //Vertex cut
    if (fabs(zvertex) > VERTEX_CUT)
        {
            return DISCARDEVENT;
        }

    hCentralityAfterVertex->Fill(centrality);
    nTotalCentralityAfterVertex += centrality;

    if (evtheader)
        {
            eventNumber = evtheader->get_EvtSequence();
        }

    //Get tracks from this event
    GetTracks(topNode, all_tracks);

    //Get clusters from this event
    GetClusters(topNode, all_clusters, zvertex);

    //Get particles from this event
    GetParticles(all_tracks, all_clusters, charged_particles, neutral_particles, all_particles);

    //Get anti-kt reconstructed jets
    //GetAntiKt(all_particles, antikt_00, antikt_jets_0, hAntikt0, hAntiktCentrality0, hAntiktErtEfficiency0, hCharged0, hNeutral0);
    GetAntiKt(all_particles, antikt_01, antikt_jets_1, hAntikt1, hAntiktCentrality1, hAntiktErtEfficiency1, hCharged1, hNeutral1);
    //GetAntiKt(all_particles, antikt_02, antikt_jets_2, hAntikt2, hAntiktCentrality2, hAntiktErtEfficiency2, hCharged2, hNeutral2);
    //GetAntiKt(all_particles, antikt_03, antikt_jets_3, hAntikt3, hAntiktCentrality3, hAntiktErtEfficiency3, hCharged3, hNeutral3);

    /////////////////////////////////////////////////////////////////////////////////////
    //Final jets
    /////////////////////////////////////////////////////////////////////////////////////
    std::vector<particles> jet_constituents;
    jet_constituents.clear();
    GetAntiKtCommon(all_particles, 0.2, 3.0, MINPT_RECO, 0.2, 0.7, final_jets, jet_constituents, 0);

    nGoodEvents[0]++;
    if(centralityBin >= 1 && centralityBin <= 4)
        {
            nGoodEvents[centralityBin]++;
        }

    for (unsigned int f = 0; f < final_jets.size(); f++)
        {
            float pT  = final_jets[f].pT;

            if(pT > 15.0)
                {
                    nJets[0]++;
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            nJets[centralityBin]++;
                        }
                }

            if(pT > 20.0)
                {
                    nJetsHighPt[0]++;
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            nJetsHighPt[centralityBin]++;
                        }
                }

            float eta = final_jets[f].eta;
            float phi = final_jets[f].phi;

            int arm = final_jets[f].arm;

            hFinalJets[0]->Fill(pT);
            hFinalJetsEtaPhi[0]->Fill(eta, phi);
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    hFinalJets[centralityBin]->Fill(pT);
                    hFinalJetsEtaPhi[centralityBin]->Fill(eta, phi);
                }
        }


    //float vNc, float vMinCf, float vMaxCf,
    //int fiducialCut, int cBin,
    makeFinalJets(all_particles,
                  3.0, 0.2, 0.7,
                  0, centralityBin,
                  hFinalJetsDefault, hFakeJetsDefault, nNoJetEvents);

    makeFinalJets(all_particles,
                  3.0, 0.2, 0.7,
                  1, centralityBin,
                  hFinalJetsEast, hFakeJetsEast, nNoJetEventsEast);

    makeFinalJets(all_particles,
                  3.0, 0.2, 0.7,
                  2, centralityBin,
                  hFinalJetsWest, hFakeJetsWest, nNoJetEventsWest);

    makeFinalJets(all_particles,
                  3.0, 0.2, 0.7,
                  3, centralityBin,
                  hFinalJetsFidTight, hFakeJetsFidTight, nNoJetEventsFidTight);

    makeFinalJets(all_particles,
                  5.0, 0.2, 0.7,
                  0, centralityBin,
                  hFinalJetsNc, hFakeJetsNc, nNoJetEventsNc);

    makeFinalJets(all_particles,
                  3.0, 0.2, 0.6,
                  0, centralityBin,
                  hFinalJetsCf, hFakeJetsCf, nNoJetEventsCf);

    makeFinalJets(all_particles,
                  5.0, 0.2, 0.6,
                  0, centralityBin,
                  hFinalJetsNcCf, hFakeJetsNcCf, nNoJetEventsNcCf);


    std::vector<particles> particles_tight;
    particles_tight.clear();
    for(int i = 0; i < all_particles.size(); i++)
        {
            int charge = all_particles[i].charge;
            float pT = all_particles[i].pT;
            float energy = all_particles[i].energy;

            bool passPt = (charge != 0) && (pT > 2.0);
            bool passEnergy = (charge == 0) && (energy > 2.0);

            if(passPt || passEnergy)
                {
                    particles_tight.push_back(all_particles[i]);
                }
        }
    makeFinalJets(particles_tight,
                  3.0, 0.2, 0.7,
                  0, centralityBin,
                  hFinalJetsTrClTight, hFakeJetsTrClTight, nNoJetEventsTrClTight);


    //Energy scale study
    //1. pi0 in EMCal
    GetPi0(all_clusters, hPi0Mass, hPi0MassVsPt, hPi0Asymmetry);

    //R=0.3
    std::vector<particles> R3_jet_constituents;
    std::vector<jets> R3_final_jets;
    GetAntiKtCommon(all_particles, 0.3, 3.0, 5.0, 0.2, 0.7,
                    R3_final_jets, R3_jet_constituents, 0);

    for (unsigned int f = 0; f < R3_final_jets.size(); f++)
        {
            float pT  = R3_final_jets[f].pT;

            if(pT > 8.0)
                {
                    hR3Jets[0]->Fill(pT);
		    if(R3_final_jets[f].arm == 0){
			hR3JetsEast->Fill(pT);
		    }else{
			hR3JetsWest->Fill(pT);
		    }
                }

            hR3AllJets[0]->Fill(pT);
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    if(pT > 8.0)
                        {
                            hR3Jets[centralityBin]->Fill(pT);
                        }

                    hR3AllJets[centralityBin]->Fill(pT);
                }


            //For ERT trigger efficiency
            bool ertTriggerFired = false;
            for (unsigned int iconst = 0; iconst < R3_jet_constituents.size(); iconst++)
                {
                    //Incase of two or more jets in an event
                    if(R3_final_jets[f].pT == R3_jet_constituents[iconst].jetPt)
                        {
                            if(R3_jet_constituents[iconst].ertTrigger == true)
                                {
                                    ertTriggerFired = true;
                                }
                        }
                }

            if(ertTriggerFired)
                {
                    float pT  = R3_final_jets[f].pT;

                    hR3ERTJets[0]->Fill(pT);
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hR3ERTJets[centralityBin]->Fill(pT);
                        }
                }
        }

    std::vector<particles> R3_shuffled_particles;
    std::vector<jets> R3_fake_jets;
    std::vector<particles> R3_fake_jet_constituents;
    if(R3_final_jets.size() == 0)
        {
            GetFakeJets(all_particles, 0.3, 3.0, 8.0, 0.2, 0.7,
                        R3_shuffled_particles, R3_fake_jets, R3_fake_jet_constituents,
                        0);

            nNoR3JetEvents[0]++;
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    nNoR3JetEvents[centralityBin]++;
                }

            for (unsigned int f = 0; f < R3_fake_jets.size(); f++)
                {
                    bool ertTriggerFired = false;
                    for (unsigned int iconst = 0; iconst < R3_fake_jet_constituents.size(); iconst++)
                        {
                            //Incase of two or more jets in an event
                            if(R3_fake_jets[f].pT == R3_fake_jet_constituents[iconst].jetPt)
                                {
                                    if(R3_fake_jet_constituents[iconst].ertTrigger == true)
                                        {
                                            ertTriggerFired = true;
                                        }
                                }
                        }
                    //For ERT triggered dataset
                    if(!isMB && !ertTriggerFired)
                        {
                            continue;
                        }

                    float pT  = R3_fake_jets[f].pT;

                    hR3FakeJets[0]->Fill(pT);
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hR3FakeJets[centralityBin]->Fill(pT);
                        }
                }
        }

    //UE study
    if(final_jets.size() == 0)
        {
            TH2F *hBackgroundEtaPhi = new TH2F("hBackgroundEtaPhi", "", 2, -0.35, 0.35, 12, -0.5, 3.7);
            for(int i = 0; i < all_particles.size(); i++)
                {
                    float eta = all_particles[i].eta;
                    float phi = all_particles[i].phi;
                    float pT = all_particles[i].pT;

                    hBackgroundEtaPhi->Fill(eta, phi, pT);
                }

            for(int ix = 1; ix <= hBackgroundEtaPhi->GetXaxis()->GetNbins(); ix++)
                {
                    for(int iy = 1; iy <= hBackgroundEtaPhi->GetYaxis()->GetNbins(); iy++)
                        {
                            if((iy == 5) || (iy == 6) || (iy == 7) || (iy == 8))
                                {
                                    continue;
                                }

                            float binContent = hBackgroundEtaPhi->GetBinContent(ix, iy);
                            hBackground[0]->Fill(binContent);
                            if(centralityBin >= 1 && centralityBin <= 4)
                                {
                                    hBackground[centralityBin]->Fill(binContent);
                                }
                        }
                }
            delete hBackgroundEtaPhi;
        }

    //Fill Trees
    FillTrees(writeTree);

    //Fill Histograms
    FillHistograms();

    // any other return code might lead to aborting the event or analysis
    return 0;
}

int JetAnalyzer::ResetEvent(PHCompositeNode *topNode)
{
    //Clear the tracks and clusters list from previous event
    all_tracks.clear();
    all_clusters.clear();

    //Clear the particle list from previous event
    charged_particles.clear();
    neutral_particles.clear();
    all_particles.clear();

    //Clear anti-kt jet list from previous event
    antikt_jets_0.clear();
    antikt_jets_1.clear();
    antikt_jets_2.clear();
    antikt_jets_3.clear();

    final_jets.clear();

    return EVENT_OK;
}

void JetAnalyzer::MyStatistics()
{
    std::ios_base::fmtflags originalFlags = std::cout.flags();

    cout.precision(3);
    cout << "++++++++++  Statistics: Tracks   +++++++++++++++++" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Number of total tracks for this run:                   " << setw(5) << nAllTracks << endl;
    cout << "Number of total tracks that pass pT:                   " << setw(5)
         << nPassPt  << setw(5) << " (" << ((float)nPassPt / nAllTracks) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass quality:              " << setw(5)
         << nPassQuality  << setw(5) << " (" << ((float)nPassQuality / nAllTracks) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass EMCal Matching:       " << setw(5)
         << nPassEMCMatching << setw(5) << " (" << ((float)nPassEMCMatching / nAllTracks) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass PC3 Matching:         " << setw(5)
         << nPassPC3Matching << setw(5) << " (" << ((float)nPassPC3Matching / nAllTracks) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass EMCal || PC3 Matching:" << setw(5)
         << nPassMatching << setw(5) << " (" << ((float)nPassMatching / nAllTracks) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass first set of cuts:    " << setw(5)
         << nPassFirst << setw(5) << " (" << ((float)nPassFirst / nAllTracks) * 100 << "%)" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    cout << "Number of total tracks that pass ghost cut:            " << setw(5)
         << nPassDcPcGhost << setw(5) << " (" << ((float)nPassDcPcGhost / nPassFirst) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass DC/Pc conversion cut: " << setw(5)
         << nPassDcPcConversion << setw(5) << " (" << ((float)nPassDcPcConversion / nPassFirst) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass Ghost and DC/PC cuts: " << setw(5)
         << nPassPairCut << setw(5) << " (" << ((float)nPassPairCut / nPassFirst) * 100 << "%)" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    cout << "Number of total tracks that pass conversion in edge:   " << setw(5)
         << nPassConversionEdge << setw(5) << " (" << ((float)nPassConversionEdge / nPassFirst) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass electron conversion:  " << setw(5)
         << nPassConversionElectron << setw(5) << " (" << ((float)nPassConversionElectron / nPassFirst) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass conversion ecore:     " << setw(5)
         << nPassConversionEcore << setw(5) << " (" << ((float)nPassConversionEcore / nPassFirst) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass conversion cut:       " << setw(5)
         << nPassConversion << setw(5) << " (" << ((float)nPassConversion / nPassFirst) * 100 << "%)" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    cout << "Number of total charged particles for this run:        " << setw(5)
         << nChargedParticles << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    cout << "++++++++++  Statistics: Clusters   +++++++++++++++++" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Number of total clusters for this run:                 " << setw(5)
         << nAllClusters << endl;
    cout << "Number of total tracks that pass isValid:              " << setw(5)
         << nPassIsValid  << setw(5) << " (" << ((float)nPassIsValid / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass notHot:               " << setw(5)
         << nPassNotHot  << setw(5) << " (" << ((float)nPassNotHot / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass notDead:              " << setw(5)
         << nPassNotDead  << setw(5) << " (" << ((float)nPassNotDead / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass notUnclaib:           " << setw(5)
         << nPassNotUncalib  << setw(5) << " (" << ((float)nPassNotUncalib / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass notBad:               " << setw(5)
         << nPassNotBad  << setw(5) << " (" << ((float)nPassNotBad / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass energy:               " << setw(5)
         << nPassEnergy  << setw(5) << " (" << ((float)nPassEnergy / nAllClusters) * 100 << "%)" << endl;
    cout << "Number of total tracks that pass everything:           " << setw(5)
         << nPassEverything  << setw(5) << " (" << ((float)nPassEverything / nAllClusters) * 100 << "%)" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Number of total neutral particles for this run:        " << setw(5)
         << nNeutralParticles << endl << endl;
    std::cout.flags(originalFlags);
}

int JetAnalyzer::EndRun(const int runNumber)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    if (verbo)
        {
            std::ios_base::fmtflags originalFlags = std::cout.flags();

            cout.precision(3);
            cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
            cout << "Run number:                                            " << runNumber << endl;
            cout << "Number of total events for this run:                   " << nTotalEvents[0] << endl;
            cout << "Number of total events for centrality 1:               " << nTotalEvents[1] << endl;
            cout << "Number of total events for centrality 2:               " << nTotalEvents[2] << endl;
            cout << "Number of total events for centrality 3:               " << nTotalEvents[3] << endl;
            cout << "Number of total events for centrality 4:               " << nTotalEvents[4] << endl << endl;
            std::cout.flags(originalFlags);

            MyStatistics();

            cout.precision(3);
            cout << "++++++++++  Statistics: Jets   +++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.2) for this run:  " << setw(5) << nJets[0] << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.2) Centrality 1:  " << setw(5) << nJets[1] << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.2) Centrality 2:  " << setw(5) << nJets[2] << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.2) Centrality 3:  " << setw(5) << nJets[3] << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.2) Centrality 4:  " << setw(5) << nJets[4] << endl << endl;
            std::cout.flags(originalFlags);
        }

    unsigned int scaledownMB = 0;
    unsigned int scaledownERT = 0;

    if(!isCuAu && isMB)
        {
            scaledownMB = triggerHelper->getLevel1Scaledown("BBCLL1(>0 tubes) narrowvtx");
            scaledownERT = triggerHelper->getLevel1Scaledown("ERT_4x4c&BBCLL1(narrow)");
        }
    nEffectiveEvents = (float)nGoodEvents[0] * (1 + scaledownMB) / (1 + scaledownERT);
    hEffectiveEvents->SetBinContent(1, nEffectiveEvents);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Good events: " << nGoodEvents[0] << endl;
    cout << "MB scaledown is: " << scaledownMB << endl;
    cout << "ERT scaledown is: " << scaledownERT << endl;
    printf("%.10E\n", nEffectiveEvents);

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    hTotalCentralityVsRun->SetBinContent(runNumber - IRUN_NUMBER + 1, nTotalCentrality);
    hCentralityEventsVsRun->SetBinContent(runNumber - IRUN_NUMBER + 1, nCentralityEvents);
    hTotalCentralityAfterVertexVsRun->SetBinContent(runNumber - IRUN_NUMBER + 1, nTotalCentralityAfterVertex);

    hEffectiveEventsVsRun->SetBinContent(runNumber - IRUN_NUMBER + 1, nEffectiveEvents);

    for(unsigned int i = 0; i < 5; i++)
        {
            hTotalVertexVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nTotalZvertex[i]);
            hVertexEventsVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nVertexEvents[i]);

            hTotalEventsVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nTotalEvents[i]);
            hGoodEventsVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nGoodEvents[i]);

            hJetYieldVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nJets[i]);
            hJetYieldHighPtVsRun[i]->SetBinContent(runNumber - IRUN_NUMBER + 1, nJetsHighPt[i]);
        }

    return EVENT_OK;
}

int JetAnalyzer::End(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  End called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    //Bins 1 to 5: Total Events
    hEvents->SetBinContent(1, nTotalEvents[0]);
    hEvents->SetBinContent(2, nTotalEvents[1]);
    hEvents->SetBinContent(3, nTotalEvents[2]);
    hEvents->SetBinContent(4, nTotalEvents[3]);
    hEvents->SetBinContent(5, nTotalEvents[4]);

    //Bins 6 to 10: Good Events
    hEvents->SetBinContent(6, nGoodEvents[0]);
    hEvents->SetBinContent(7, nGoodEvents[1]);
    hEvents->SetBinContent(8, nGoodEvents[2]);
    hEvents->SetBinContent(9, nGoodEvents[3]);
    hEvents->SetBinContent(10, nGoodEvents[4]);

    //For fake jets
    hEvents->SetBinContent(11, nNoJetEvents[0]);
    hEvents->SetBinContent(12, nNoJetEvents[1]);
    hEvents->SetBinContent(13, nNoJetEvents[2]);
    hEvents->SetBinContent(14, nNoJetEvents[3]);
    hEvents->SetBinContent(15, nNoJetEvents[4]);

    hEvents->SetBinContent(16, nNoJetEventsEast[0]);
    hEvents->SetBinContent(17, nNoJetEventsEast[1]);
    hEvents->SetBinContent(18, nNoJetEventsEast[2]);
    hEvents->SetBinContent(19, nNoJetEventsEast[3]);
    hEvents->SetBinContent(20, nNoJetEventsEast[4]);

    hEvents->SetBinContent(21, nNoJetEventsWest[0]);
    hEvents->SetBinContent(22, nNoJetEventsWest[1]);
    hEvents->SetBinContent(23, nNoJetEventsWest[2]);
    hEvents->SetBinContent(24, nNoJetEventsWest[3]);
    hEvents->SetBinContent(25, nNoJetEventsWest[4]);

    hEvents->SetBinContent(26, nNoJetEventsFidTight[0]);
    hEvents->SetBinContent(27, nNoJetEventsFidTight[1]);
    hEvents->SetBinContent(28, nNoJetEventsFidTight[2]);
    hEvents->SetBinContent(29, nNoJetEventsFidTight[3]);
    hEvents->SetBinContent(30, nNoJetEventsFidTight[4]);

    hEvents->SetBinContent(31, nNoJetEventsNc[0]);
    hEvents->SetBinContent(32, nNoJetEventsNc[1]);
    hEvents->SetBinContent(33, nNoJetEventsNc[2]);
    hEvents->SetBinContent(34, nNoJetEventsNc[3]);
    hEvents->SetBinContent(35, nNoJetEventsNc[4]);

    hEvents->SetBinContent(36, nNoJetEventsCf[0]);
    hEvents->SetBinContent(37, nNoJetEventsCf[1]);
    hEvents->SetBinContent(38, nNoJetEventsCf[2]);
    hEvents->SetBinContent(39, nNoJetEventsCf[3]);
    hEvents->SetBinContent(40, nNoJetEventsCf[4]);

    hEvents->SetBinContent(41, nNoJetEventsNcCf[0]);
    hEvents->SetBinContent(42, nNoJetEventsNcCf[1]);
    hEvents->SetBinContent(43, nNoJetEventsNcCf[2]);
    hEvents->SetBinContent(44, nNoJetEventsNcCf[3]);
    hEvents->SetBinContent(45, nNoJetEventsNcCf[4]);

    hEvents->SetBinContent(46, nNoJetEventsTrClTight[0]);
    hEvents->SetBinContent(47, nNoJetEventsTrClTight[1]);
    hEvents->SetBinContent(48, nNoJetEventsTrClTight[2]);
    hEvents->SetBinContent(49, nNoJetEventsTrClTight[3]);
    hEvents->SetBinContent(50, nNoJetEventsTrClTight[4]);

    //R=0.3 fakes
    hEvents->SetBinContent(60, nNoR3JetEvents[0]);
    hEvents->SetBinContent(61, nNoR3JetEvents[1]);
    hEvents->SetBinContent(62, nNoR3JetEvents[2]);
    hEvents->SetBinContent(63, nNoR3JetEvents[3]);
    hEvents->SetBinContent(64, nNoR3JetEvents[4]);

    outfile->Write();
    outfile->Close();

    delete antikt_00;
    delete antikt_01;
    delete antikt_02;
    delete antikt_03;

    delete triggerHelper;

    return EVENT_OK;
}

void JetAnalyzer::SetData(bool what = true)
{
    isData = what;
}

void JetAnalyzer::SetCuAu(bool what = true)
{
    cout << endl << endl;
    cout << "***********************************************************************" << endl;
    cout << "***********************************************************************" << endl;
    if(what)
        {
            cout << endl;
            cout << "Running over Cu+Au data setup" << endl << endl;
        }
    else
        {
            cout << endl;
            cout << "Running over p+p data setup" << endl << endl;
        }

    isCuAu = what;
    cout << "***********************************************************************" << endl << endl;
}

void JetAnalyzer::SetMB(bool what = true)
{
    if(what)
        {
            cout << "Running over MB data" << endl << endl;

            if(isCuAu)
                {
                    cout << "Looking for Min Bias: BBCLL1(>1 tubes) narrowvtx, bit mask: " << MB_CUAU_NARROWVTX << endl << endl;
                }
            else
                {
                    cout << "Looking for Min Bias: BBCLL1(>0 tubes) narrowvtx, bit mask: " << MB_PP_NARROWVTX << endl << endl;
                }
        }
    else
        {
            cout << "Running over ERT data" << endl << endl;
            cout << "Looking for ERT 4x4c: ERT_4x4c&BBCLL1(narrow), bit mask: " << ERT_BIT_MASK_C << endl << endl;
        }

    isMB = what;
    cout << "***********************************************************************" << endl;
    cout << "***********************************************************************" << endl << endl;
}

void JetAnalyzer::MyInit()
{
    if(isCuAu)
        {
            //Load and read the warnmap for Cu+Au
            EmcMap::ReadWarnMap("warnmapCuAu.txt");

            //For PC3 mathing
            pc3_init_fit_pars_II();

            //For EMC matching
            LoadEmcMatchingParameters();
        }
    else
        {
            //Load and read the warnmap for P+P
            EmcMap::ReadWarnMap("warnmapPP.txt");

            //For p+p Quality map
            TrackQualityPP::InitTrackQuality();
        }
}

void JetAnalyzer::GetTracks(PHCompositeNode *topNode, std::vector<tracks>& track_list)
{
    PHCentralTrack  *phcentral = getClass<PHCentralTrack>(topNode, "PHCentralTrack");
    if (!phcentral)
        {
            cout << "No PHCentral!  No sense continuing" << endl;
            exit(1);
        }

    for (unsigned int i = 0; i < phcentral->get_npart(); i++)
        {
            float mom           = phcentral->get_mom(i); // Magnitude of the momentum.
            float theta         = phcentral->get_the0(i); //The track's theta direction at the vertex
            float pT            = mom * sin(theta);
            float eta           = (float)(-log(tan(0.5 * theta)));
            float phi           = phcentral->get_phi0(i); // The track's phi direction at the vertex
            float phiDC         = phcentral->get_phi(i); // The phi coordinate at which the track crosses the drift chamber reference radius
            float zedDC         = phcentral->get_zed(i); //The Z coordinate at which the track crosses PC1
            float alpha         = phcentral->get_alpha(i); // This is the inclination of the track w.r.t. an infinite momentum track
            float energy        = phcentral->get_ecore(i); // EMC "shower core" energy

            int charge          = phcentral->get_charge(i);
            int quality         = phcentral->get_quality(i); // Quality of the Drift Chamber Tracks
            int n0              = phcentral->get_n0(i); //The number of phototubes that fired in the normally sized ring area
            int arm             = phcentral->get_dcarm(i); //Arm containing the track (East=0, West=1)
            int sector          = phcentral->get_sect(i); //EMC sector of the associated cluster
            int emcid           = phcentral->get_emcid(i); // Index of the emc cluster used

            float pc3dphi       = phcentral->get_pc3dphi(i);
            float pc3sdphi      = phcentral->get_pc3sdphi(i);
            float pc3dz         = phcentral->get_pc3dz(i);
            float pc3sdz        = phcentral->get_pc3sdz(i);
            float emcdphi       = phcentral->get_emcdphi(i);
            float emcsdphi      = phcentral->get_emcsdphi(i);
            float emcdz         = phcentral->get_emcdz(i);
            float emcsdz        = phcentral->get_emcsdz(i);

            float beta          = phcentral->get_beta(i);

            nAllTracks++;

            if(isCuAu)
                {
                    float newAlpha  = getNewAlpha(alpha, phiDC, arm);// phi is phiDC not phi0 !!!!!!!!!!!
                    float deltaAlpha = newAlpha - alpha;

                    pT = pT * fabs(alpha / newAlpha);
                    phi = 2.0195 * deltaAlpha + phi;

                    //change the pt and momentum
                    if (arm == 0)
                        {
                            pT = pT * MomfactorEast;
                        }
                    if (arm == 1)
                        {
                            pT = pT * MomfactorWest;
                        }
                    mom = pT / sin(theta);

                    alpha = newAlpha;
                }

            float px = pT * cos(phi);
            float py = pT * sin(phi);
            float pz = mom * cos(theta);
            float eT = energy * sin(theta);

            int armsect = -99;
            if (arm == 1)
                {
                    armsect = sector;
                }
            else
                {
                    armsect = 4 + sector;
                }

            //Maximum and minimum momentum cut- higher momentum likely to come from conversions. Modified Oct 2, 2013- just minimum cut
            bool passPt = (pT > TRACK_MIN_PT_CUT);
            if (passPt)
                {
                    nPassPt++;
                }

            //For Quality Mask
            bool x1Used = ((quality & X1_USED) == X1_USED);
            bool x2Used = ((quality & X2_USED) == X2_USED);
            bool uvUsed = ((quality & UV_FOUND_UNIQUE) == UV_FOUND_UNIQUE);

            float board = -9999.99;
            bool inBrokenX1 = false;
            bool inBrokenX2 = false;
            bool inBrokenUV = false;
            bool passQuality = false;

            if(isCuAu)
                {
                    board = TrackQualityCuAu::getBoard(phiDC, arm);
                    inBrokenX1 = TrackQualityCuAu::inBrokenX1(phiDC, alpha, zedDC, arm);
                    inBrokenX2 = TrackQualityCuAu::inBrokenX2(phiDC, alpha, zedDC, arm);
                    inBrokenUV = TrackQualityCuAu::inBrokenUV(phiDC, alpha, zedDC, arm);

                    //Quality cut
                    passQuality = TrackQualityCuAu::passQualityMask(quality, phiDC, alpha, zedDC, arm);
                }
            else
                {
                    board = TrackQualityPP::getBoard(phiDC, arm);
                    inBrokenX1 = TrackQualityPP::inBrokenX1(phiDC, alpha, zedDC, arm);
                    inBrokenX2 = TrackQualityPP::inBrokenX2(phiDC, alpha, zedDC, arm);
                    inBrokenUV = TrackQualityPP::inBrokenUV(phiDC, alpha, zedDC, arm);

                    //Quality cut
                    passQuality = TrackQualityPP::passQualityMask(quality, phiDC, alpha, zedDC, arm);
                }

            if (passQuality)
                {
                    nPassQuality++;
                }

            //To be consistant between p+p and Cu+Au dataset + to match with simulation
            bool passZedDC = fabs(zedDC) < 75 && fabs(zedDC) > 3;

            bool passMost = passPt && armsect >= 0 && armsect <= 7;
            bool passDC = passMost && passQuality && passZedDC;

            if(isCuAu)
                {
                    //For PC3 matching
                    int ized = -9999;
                    if     (zedDC > -70 && zedDC <= -57)
                        {
                            ized = 0;
                        }
                    else if (zedDC > -57 && zedDC <= -44)
                        {
                            ized = 1;
                        }
                    else if (zedDC > -44 && zedDC <= -31)
                        {
                            ized = 2;
                        }
                    else if (zedDC > -31 && zedDC <= -18)
                        {
                            ized = 3;
                        }
                    else if (zedDC > -18 && zedDC <= -5)
                        {
                            ized = 4;
                        }
                    else if (zedDC >= 5  && zedDC <  18)
                        {
                            ized = 5;
                        }
                    else if (zedDC >= 18 && zedDC <  31)
                        {
                            ized = 6;
                        }
                    else if (zedDC >= 31 && zedDC <  44)
                        {
                            ized = 7;
                        }
                    else if (zedDC >= 44 && zedDC <  57)
                        {
                            ized = 8;
                        }
                    else if (zedDC >= 57 && zedDC <  70)
                        {
                            ized = 9;
                        }

                    if (!(ized < 0 || ized > 9) && passMost)
                        {
                            int NMUL = 10;
                            int icent = (int) ( NMUL * ((centrality - 0.001) / 100.) );
                            if (!(icent < 0 || icent > NMUL))
                                {
                                    pc3sdphi = pc3_sdphi_func_II(charge, arm, icent, ized, pT, pc3dphi);
                                    pc3sdz   = pc3_sdz_func_II(charge, arm, icent, ized, pT, pc3dz);
                                }
                        }

                    //For EMC matching
                    if(fabs(zedDC) < 75)
                        {
                            int zbin = (int)(NZED * (zedDC + 75.0) / 150.0);
                            int chargeBin = 1;
                            if (charge == -1)
                                {
                                    chargeBin = 0;
                                }

                            int cbin = -9999;
                            if (centrality > 0 && centrality <= 10)
                                {
                                    cbin = 0;
                                }
                            if (centrality > 10 && centrality <= 20)
                                {
                                    cbin = 1;
                                }
                            if (centrality > 20 && centrality <= 30)
                                {
                                    cbin = 2;
                                }
                            if (centrality > 30 && centrality <= 40)
                                {
                                    cbin = 3;
                                }
                            if (centrality > 40 && centrality <= 60)
                                {
                                    cbin = 4;
                                }
                            if (centrality > 60 && centrality <= 100)
                                {
                                    cbin = 5;
                                }

                            if ((zbin >= 0 && zbin <= 9) && (cbin >= 0 && cbin <= 5) && passMost)
                                {
                                    float emcsdphi_initial = CalculateInitialEmcsdPhi(chargeBin, armsect, zbin, cbin, emcdphi, pT);
                                    emcsdphi = CalculateFinalEmcsdPhi(chargeBin, armsect, zbin, cbin, emcsdphi_initial, pT);

                                    float emcdz_corrected = CalculateCorrectedEmcdZ(beta, pT, emcdz, armsect);
                                    float emcsdz_initial = CalculateInitialEmcsdZ(chargeBin, armsect, zbin, cbin, emcdz_corrected, pT);
                                    emcsdz = CalculateFinalEmcsdZ(chargeBin, armsect, zbin, cbin, emcsdz_initial, pT);
                                }
                        }
                }

            //Matching cut
            bool emcMatching = sqrt((emcsdphi * emcsdphi) + (emcsdz * emcsdz)) < 3.0;
            bool pc3Matching = sqrt((pc3sdphi * pc3sdphi) + (pc3sdz * pc3sdz)) < 3.0;

            if(emcMatching)
                {
                    nPassEMCMatching++;
                }

            if(pc3Matching)
                {
                    nPassPC3Matching++;
                }

            bool passMatching = pc3Matching || emcMatching;

            if (passMatching)
                {
                    nPassMatching++;
                }

            bool passFirst = passPt && passQuality && passMatching && passZedDC;

            if (passFirst)
                {
                    nPassFirst++;
                }

            tracks temp;
            temp.mom            = mom;
            temp.theta          = theta;
            temp.pT             = pT;
            temp.eta            = eta;
            temp.phi            = phi;
            temp.phiDC          = phiDC;
            temp.zedDC          = zedDC;
            temp.alpha          = alpha;
            temp.energy         = energy;
            temp.board          = board;

            temp.charge         = charge;
            temp.quality        = quality;
            temp.n0             = n0;
            temp.arm            = arm;
            temp.armsect        = armsect;
            temp.emcid          = emcid;

            temp.px             = px;
            temp.py             = py;
            temp.pz             = pz;
            temp.eT             = eT;

            temp.pc3dphi        = pc3dphi;
            temp.pc3dz          = pc3dz;
            temp.pc3sdphi       = pc3sdphi;
            temp.pc3sdz         = pc3sdz;
            temp.emcdphi        = emcdphi;
            temp.emcdz          = emcdz;
            temp.emcsdphi       = emcsdphi;
            temp.emcsdz         = emcsdz;

            temp.inBrokenX1     = inBrokenX1;
            temp.inBrokenX2     = inBrokenX2;
            temp.inBrokenUV     = inBrokenUV;
            temp.x1Used         = x1Used;
            temp.x2Used         = x2Used;
            temp.uvUsed         = uvUsed;
            temp.passDC         = passDC;
            temp.passQuality    = passQuality;
            temp.passMatching   = passMatching;
            temp.passFirst      = passFirst;

            track_list.push_back(temp);
        }
}


void JetAnalyzer::GetClusters(PHCompositeNode *topNode, std::vector<clusters>& cluster_list, float vertex, bool perfectEMCal)
{
    //The number of towers is 48(y) x 96(z) x 2(PbGl) + 36(y) x 72(z) x 6(PbSc) = 9216+15552 = 24768.
    //Tower ID is a unique ID of all towers in eight sectors (0~24767).
    emcClusterContainer *emcclustercontainer = getClass<emcClusterContainer>(topNode, "emcClusterContainer");
    if (!emcclustercontainer)
        {
            cout << "No emcClusterContainer!  No sense continuing" << endl;
            exit(1);
        }

    if(isData)
        {
            //For ERT trigger
            ertOut = getClass<ErtOut>(topNode, "ErtOut");
            if (!ertOut)
                {
                    cout << "No ertOut!  No sense continuing" << endl;
                    exit(1);
                }
        }

    int Nclus = emcclustercontainer->size();
    for (int iclus = 0; iclus < Nclus; iclus++)
        {
            nAllClusters++;

            emcClusterContent* clus = emcclustercontainer->getCluster(iclus);
            int arm         = clus->arm(); //In EMCal convention, West Arm is 0 and East Arm is 1, and thus armsector 0...7 are W0...W3 E0...E3
            int sector      = clus->sector();
            int emcid       = clus->id();
            float energy    = clus->ecore();
            float prob      = clus->prob_photon();
            float chi2      = clus->chi2();

            //need to construct cluster theta and pT
            float x         = clus->x();
            float y         = clus->y();
            float z         = clus->z() - vertex;
            float theta     = clus->theta();//acos(z / sqrt(x * x + y * y + z * z));

            TVector3 v3;
            TLorentzVector v4;
            v3.SetXYZ(x, y, z);
            v3 = energy * v3.Unit();
            v4.SetT(energy);
            v4.SetVect(v3);

            float pT        = v4.Pt();
            float phi       = phiReduce(v4.Phi());
            float px        = v4.Px();
            float py        = v4.Py();
            float pz        = v4.Pz();
            float eta       = v4.Eta();
            float eT        = v4.Et();

            int yTowerPos   = clus->iypos();
            int zTowerPos   = clus->izpos();

            int armsect     = (arm * 4) + sector; //Arm for emcal is different

            int towerId         = EmcMap::GetTowerID(armsect, yTowerPos, zTowerPos);
            bool passIsValid    = EmcMap::IsValidYZ(armsect, yTowerPos, zTowerPos);
            bool passNotHot     = !EmcMap::IsHot(armsect, yTowerPos, zTowerPos);
            bool passNotDead    = !EmcMap::IsDead(armsect, yTowerPos, zTowerPos);
            bool passNotUncalib = !EmcMap::IsUncalib(armsect, yTowerPos, zTowerPos);
            bool passNotBad     = !EmcMap::IsBad(armsect, yTowerPos, zTowerPos);

            //ERT information
            int ertSM = (arm == 1 && sector < 2)  ?  (yTowerPos / 12) * 8 + zTowerPos / 12  :  (yTowerPos / 12) * 6 + zTowerPos / 12;

            int ert4x4Bit = 0;
            if(isData)
                {
                    //ertBit: 0=4x4a, 1=4x4b, 2=4x4c
                    ert4x4Bit = ertOut->get_ERTbit(2, arm, sector, ertSM);
                }

            if (passIsValid)
                {
                    nPassIsValid++;
                }
            if (passNotHot)
                {
                    nPassNotHot++;
                }
            if (passNotDead)
                {
                    nPassNotDead++;
                }
            if (passNotUncalib)
                {
                    nPassNotUncalib++;
                }
            if (passNotBad)
                {
                    nPassNotBad++;
                }

            bool passEnergy = energy > CLUSTER_MIN_ENERGY_CUT;
            if (passEnergy)
                {
                    nPassEnergy++;
                }

            bool passEverything = passIsValid && passNotBad && passEnergy;

            //Needed for simulation study
            if(perfectEMCal)
                {
                    passEverything = passEnergy;
                }

            if (passEverything)
                {
                    nPassEverything++;
                }

            //Convert EMC arm to DC convention- need this for simulation
            if(arm == 0)
                {
                    arm = 1;
                }
            else
                {
                    arm = 0;
                }
            clusters temp;
            temp.arm             = arm;
            temp.sector          = sector;
            temp.armsect         = armsect;
            temp.emcid           = emcid;
            temp.yTowerPos       = yTowerPos;
            temp.zTowerPos       = zTowerPos;
            temp.towerId         = towerId;
            temp.ert4x4Bit       = ert4x4Bit;

            temp.energy          = energy;
            temp.theta           = theta;
            temp.pT              = pT;
            temp.eT              = eT;
            temp.phi             = phi;
            temp.px              = px;
            temp.py              = py;
            temp.pz              = pz;
            temp.eta             = eta;
            temp.chi2            = chi2;
            temp.prob            = prob;

            temp.passIsValid     = passIsValid;
            temp.passNotHot      = passNotHot;
            temp.passNotDead     = passNotDead;
            temp.passNotUncalib  = passNotUncalib;
            temp.passNotBad      = passNotBad;
            temp.passEnergy      = passEnergy;
            temp.passEverything  = passEverything;

            cluster_list.push_back(temp);
        }
}

void JetAnalyzer::GetParticles(std::vector<tracks> track_list,
                               std::vector<clusters> cluster_list,
                               std::vector<tracks>& charged_particle_list,
                               std::vector<clusters>& neutral_particle_list,
                               std::vector<particles>& particle_list)
{
    //Ghost and conversion cut in DC/PC1
    unsigned int indexPair;
    if(track_list.size() != 0)
        {
            indexPair = track_list.size();
        }
    else
        {
            indexPair = 1;
        }

    bool passGhost[indexPair];
    fill(passGhost, passGhost + indexPair / sizeof(bool), true);
    bool passDcConversion[indexPair];
    fill(passDcConversion, passDcConversion + indexPair / sizeof(bool), true);

    for (unsigned int t = 0; t < track_list.size(); t++)
        {

            if (!track_list[t].passFirst)
                {
                    passGhost[t] = false;
                    passDcConversion[t] = false;
                    continue;
                }
            for (unsigned int tt = t + 1; tt < track_list.size(); tt++)
                {
                    if (!track_list[tt].passFirst)
                        {
                            continue;
                        }

                    if (track_list[t].arm != track_list[tt].arm)
                        {
                            continue;
                        }
                    float dPhi = track_list[t].phiDC - track_list[tt].phiDC;
                    float dZed = track_list[t].zedDC - track_list[tt].zedDC;

                    //Ghost cut: Ghost pairs have same charge
                    if (track_list[t].charge == track_list[tt].charge && fabs(dPhi) < 0.024 && fabs(dZed) < 0.105)
                        {
                            float asymmetry = fabs(track_list[t].pT - track_list[tt].pT) / (track_list[t].pT + track_list[tt].pT);

                            //If assymetry in pT is <0.3 reject only one track- else reject both
                            if (asymmetry < 0.3)
                                {
                                    passGhost[tt] = false;
                                }
                            else
                                {
                                    passGhost[t] = false;
                                    passGhost[tt] = false;
                                }
                        }

                    //Converion in DC/PC1: Conversion pairs have different charge
                    if (track_list[t].charge != track_list[tt].charge && fabs(dPhi) < 0.07 && fabs(dZed) < 0.105)
                        {
                            passDcConversion[t] = false;
                            passDcConversion[tt] = false;
                        }
                }
        }

    //We obtain differnt particles from tracks and clusters here
    for (unsigned int t = 0; t < track_list.size(); t++)
        {
            float mom            = track_list[t].mom;
            float theta          = track_list[t].theta;
            float pT             = track_list[t].pT;
            float eta            = track_list[t].eta;
            float phi            = track_list[t].phi;
            float phiDC          = track_list[t].phiDC;
            float zedDC          = track_list[t].zedDC;
            float alpha          = track_list[t].alpha;
            float energy         = track_list[t].energy;
            float board          = track_list[t].board;

            int charge           = track_list[t].charge;
            int quality          = track_list[t].quality;
            int n0               = track_list[t].n0;
            int arm              = track_list[t].arm;
            int armsect          = track_list[t].armsect;
            int emcid            = track_list[t].emcid;

            float px             = track_list[t].px;
            float py             = track_list[t].py;
            float pz             = track_list[t].pz;
            float eT             = track_list[t].eT;

            float pc3dphi        = track_list[t].pc3dphi;
            float pc3dz          = track_list[t].pc3dz;
            float pc3sdphi       = track_list[t].pc3sdphi;
            float pc3sdz         = track_list[t].pc3sdz;
            float emcdphi        = track_list[t].emcdphi;
            float emcdz          = track_list[t].emcdz;
            float emcsdphi       = track_list[t].emcsdphi;
            float emcsdz         = track_list[t].emcsdz;

            bool inBrokenX1     = track_list[t].inBrokenX1;
            bool inBrokenX2     = track_list[t].inBrokenX2;
            bool inBrokenUV     = track_list[t].inBrokenUV;
            bool x1Used         = track_list[t].x1Used;
            bool x2Used         = track_list[t].x2Used;
            bool uvUsed         = track_list[t].uvUsed;
            bool passDC         = track_list[t].passDC;
            bool passQuality    = track_list[t].passQuality;
            bool passMatching   = track_list[t].passMatching;
            bool passFirst      = track_list[t].passFirst;

            if(!passFirst)
                {
                    continue;
                }

            //Additional cuts:
            //Track should not be a ghost
            bool passDcPcGhosts = passGhost[t];
            if (passDcPcGhosts)
                {
                    nPassDcPcGhost++;
                }

            //Track should not be from conversion in DC/PC1
            bool passDcPcConversions = passDcConversion[t];
            if (passDcPcConversions)
                {
                    nPassDcPcConversion++;
                }

            bool passPair = passDcPcGhosts && passDcPcConversions;
            if (passPair)
                {
                    nPassPairCut++;
                }

            bool emcMatching = sqrt((emcsdphi * emcsdphi) + (emcsdz * emcsdz)) < 3.0;

            //Conversion cuts
            bool conversionEdge = inEdge(phi);
            if(!conversionEdge)
                {
                    nPassConversionEdge++;
                }

            bool conversionElectron = pT < 4.5 && n0 >= 2 && (energy / mom) < 0.6;
            if(!conversionElectron)
                {
                    nPassConversionElectron++;
                }

            bool conversionEcore = emcMatching && energy < 0.2;
            if(!conversionEcore)
                {
                    nPassConversionEcore++;
                }

            bool passConversions = !conversionEdge && !conversionElectron && !conversionEcore;
            if (passConversions)
                {
                    nPassConversion++;
                }

            bool passEverything = passPair && passConversions;

            if (!passEverything)
                {
                    continue;
                }

            nChargedParticles++;

            tracks temp;
            temp.mom            = mom;
            temp.theta          = theta;
            temp.pT             = pT;
            temp.eta            = eta;
            temp.phi            = phi;
            temp.phiDC          = phiDC;
            temp.zedDC          = zedDC;
            temp.alpha          = alpha;
            temp.energy         = energy;
            temp.board          = board;

            temp.charge         = charge;
            temp.quality        = quality;
            temp.n0             = n0;
            temp.arm            = arm;
            temp.armsect        = armsect;
            temp.emcid          = emcid;

            temp.px             = px;
            temp.py             = py;
            temp.pz             = pz;
            temp.eT             = eT;

            temp.pc3dphi        = pc3dphi;
            temp.pc3dz          = pc3dz;
            temp.pc3sdphi       = pc3sdphi;
            temp.pc3sdz         = pc3sdz;
            temp.emcdphi        = emcdphi;
            temp.emcdz          = emcdz;
            temp.emcsdphi       = emcsdphi;
            temp.emcsdz         = emcsdz;

            temp.inBrokenX1     = inBrokenX1;
            temp.inBrokenX2     = inBrokenX2;
            temp.inBrokenUV     = inBrokenUV;
            temp.x1Used         = x1Used;
            temp.x2Used         = x2Used;
            temp.uvUsed         = uvUsed;
            temp.passDC         = passDC;
            temp.passQuality    = passQuality;
            temp.passMatching   = passMatching;
            temp.passFirst      = passFirst;

            charged_particle_list.push_back(temp);
        }


    //Neutral
    for (unsigned int c = 0; c < cluster_list.size(); c++)
        {
            //Pass Everything
            if (!cluster_list[c].passEverything)
                {
                    continue;
                }

            //Don't double count
            bool charged = false;
            for (unsigned int t = 0; t < track_list.size(); t++)
                {
                    if(!track_list[t].passDC)
                        {
                            continue;
                        }
                    float emcsdphi = track_list[t].emcsdphi;
                    float emcsdz   = track_list[t].emcsdz;

                    bool emcMatching = sqrt((emcsdphi * emcsdphi) + (emcsdz * emcsdz)) < 3.0;

                    if ((cluster_list[c].emcid == track_list[t].emcid) && emcMatching)
                        {
                            charged = true;
                        }
                }

            if (!charged)
                {
                    nNeutralParticles++;

                    clusters temp;
                    temp.arm             = cluster_list[c].arm;
                    temp.sector          = cluster_list[c].sector;
                    temp.armsect         = cluster_list[c].armsect;
                    temp.emcid           = cluster_list[c].emcid;
                    temp.yTowerPos       = cluster_list[c].yTowerPos;
                    temp.zTowerPos       = cluster_list[c].zTowerPos;
                    temp.towerId         = cluster_list[c].towerId;
                    temp.ert4x4Bit       = cluster_list[c].ert4x4Bit;

                    temp.energy          = cluster_list[c].energy;
                    temp.theta           = cluster_list[c].theta;
                    temp.pT              = cluster_list[c].pT;
                    temp.eT              = cluster_list[c].eT;
                    temp.phi             = cluster_list[c].phi;
                    temp.px              = cluster_list[c].px;
                    temp.py              = cluster_list[c].py;
                    temp.pz              = cluster_list[c].pz;
                    temp.eta             = cluster_list[c].eta;
                    temp.chi2            = cluster_list[c].chi2;
                    temp.prob            = cluster_list[c].prob;

                    temp.passIsValid     = cluster_list[c].passIsValid;
                    temp.passNotHot      = cluster_list[c].passNotHot;
                    temp.passNotDead     = cluster_list[c].passNotDead;
                    temp.passNotUncalib  = cluster_list[c].passNotUncalib;
                    temp.passNotBad      = cluster_list[c].passNotBad;
                    temp.passEnergy      = cluster_list[c].passEnergy;
                    temp.passEverything  = cluster_list[c].passEverything;

                    neutral_particle_list.push_back(temp);
                }
        }

    unsigned int pId = 0;
    particles tempParticle;
    for(unsigned int c = 0; c < charged_particle_list.size(); c++)
        {
            tempParticle.arm        = charged_particle_list[c].arm;
            tempParticle.charge     = charged_particle_list[c].charge;

            tempParticle.energy     = charged_particle_list[c].energy;
            tempParticle.mom        = charged_particle_list[c].mom;
            tempParticle.pT         = charged_particle_list[c].pT;
            tempParticle.eT         = charged_particle_list[c].eT;
            tempParticle.px         = charged_particle_list[c].px;
            tempParticle.py         = charged_particle_list[c].py;
            tempParticle.pz         = charged_particle_list[c].pz;
            tempParticle.eta        = charged_particle_list[c].eta;
            tempParticle.phi        = charged_particle_list[c].phi;
            tempParticle.phiDC      = charged_particle_list[c].phiDC;

            tempParticle.ertTrigger = false;

            tempParticle.id         = pId;
            tempParticle.jetPt      = -999.9;

            particle_list.push_back(tempParticle);

            pId++;
        }

    for(unsigned int n = 0; n < neutral_particle_list.size(); n++)
        {
            tempParticle.arm        = neutral_particle_list[n].arm;
            tempParticle.charge     = 0;

            tempParticle.energy     = neutral_particle_list[n].energy;
            tempParticle.mom        = sqrt((neutral_particle_list[n].px * neutral_particle_list[n].px) +
                                           (neutral_particle_list[n].py * neutral_particle_list[n].py) +
                                           (neutral_particle_list[n].pz * neutral_particle_list[n].pz));
            tempParticle.pT         = neutral_particle_list[n].pT;
            tempParticle.eT         = neutral_particle_list[n].eT;
            tempParticle.px         = neutral_particle_list[n].px;
            tempParticle.py         = neutral_particle_list[n].py;
            tempParticle.pz         = neutral_particle_list[n].pz;
            tempParticle.eta        = neutral_particle_list[n].eta;
            tempParticle.phi        = neutral_particle_list[n].phi;
            tempParticle.phiDC      = neutral_particle_list[n].phi;

            tempParticle.ertTrigger = neutral_particle_list[n].ert4x4Bit;

            tempParticle.id         = pId;
            tempParticle.jetPt      = -999.9;

            particle_list.push_back(tempParticle);

            pId++;
        }
}

void JetAnalyzer::GetAntiKtCommon(std::vector<particles> particle_list,
                                  float R, float nc, float minPt, float minCf, float maxCf,
                                  std::vector<jets>& antikt_jet_list,
                                  std::vector<particles>& constituent_list,
                                  int fiducialCut)
{
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    //Generate input for jet reconstuction
    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    //Sort the particle list- charged particles first and then neurtal
    std::sort(particle_list.begin(), particle_list.end(), sortParticle());

    unsigned int indexTotal = 0;
    unsigned int indexCharged = 0;
    if(particle_list.size() != 0)
        {
            for (unsigned int h = 0; h < particle_list.size(); h++)
                {
                    indexTotal++;
                    if (particle_list[h].charge != 0)
                        {
                            indexCharged++;
                        }
                }
        }
    else
        {
            indexTotal = 1;
        }


    int particleArm[indexTotal];
    fill(particleArm, particleArm + indexTotal / sizeof(int), -1);
    int particleCharge[indexTotal];
    fill(particleCharge, particleCharge + indexTotal / sizeof(int), -1);

    float particleEnergy[indexTotal];
    fill(particleEnergy, particleEnergy + indexTotal / sizeof(float), -999.9);
    float particleMom[indexTotal];
    fill(particleMom, particleMom + indexTotal / sizeof(float), -999.9);
    float particlePt[indexTotal];
    fill(particlePt, particlePt + indexTotal / sizeof(float), -999.9);
    float particleEt[indexTotal];
    fill(particleEt, particleEt + indexTotal / sizeof(float), -999.9);
    float particlePx[indexTotal];
    fill(particlePx, particlePx + indexTotal / sizeof(float), -999.9);
    float particlePy[indexTotal];
    fill(particlePy, particlePy + indexTotal / sizeof(float), -999.9);
    float particlePz[indexTotal];
    fill(particlePz, particlePz + indexTotal / sizeof(float), -999.9);
    float particleEta[indexTotal];
    fill(particleEta, particleEta + indexTotal / sizeof(float), -999.9);
    float particlePhi[indexTotal];
    fill(particlePhi, particlePhi + indexTotal / sizeof(float), -999.9);
    float particlePhiDC[indexTotal];
    fill(particlePhiDC, particlePhiDC + indexTotal / sizeof(float), -999.9);

    bool particleErtTrigger[indexTotal];
    fill(particleErtTrigger, particleErtTrigger + indexTotal / sizeof(bool), false);

    int particleId[indexTotal];
    fill(particleId, particleId + indexTotal / sizeof(int), -1);

    int index = 0;
    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            if (particle_list[h].charge != 0)
                {
                    fastjet::PseudoJet pseudoCharged(particle_list[h].px,
                                                     particle_list[h].py,
                                                     particle_list[h].pz,
                                                     particle_list[h].mom);
                    pseudoCharged.set_user_index(index);

                    particleArm[index]        = particle_list[h].arm;
                    particleCharge[index]     = particle_list[h].charge;

                    particleEnergy[index]     = particle_list[h].energy;
                    particleMom[index]        = particle_list[h].mom;
                    particlePt[index]         = particle_list[h].pT;
                    particleEt[index]         = particle_list[h].eT;
                    particlePx[index]         = particle_list[h].px;
                    particlePy[index]         = particle_list[h].py;
                    particlePz[index]         = particle_list[h].pz;
                    particleEta[index]        = particle_list[h].eta;
                    particlePhi[index]        = particle_list[h].phi;
                    particlePhiDC[index]      = particle_list[h].phiDC;

                    particleErtTrigger[index] = particle_list[h].ertTrigger;

                    particleId[index]         = particle_list[h].id;

                    jetParticles_all.push_back(pseudoCharged);
                    index++;
                }
            else
                {
                    fastjet::PseudoJet pseudoNeutral(particle_list[h].px,
                                                     particle_list[h].py,
                                                     particle_list[h].pz,
                                                     particle_list[h].energy);
                    pseudoNeutral.set_user_index(index);

                    particleArm[index]        = particle_list[h].arm;
                    particleCharge[index]     = particle_list[h].charge;

                    particleEnergy[index]     = particle_list[h].energy;
                    particleMom[index]        = particle_list[h].mom;
                    particlePt[index]         = particle_list[h].pT;
                    particleEt[index]         = particle_list[h].eT;
                    particlePx[index]         = particle_list[h].px;
                    particlePy[index]         = particle_list[h].py;
                    particlePz[index]         = particle_list[h].pz;
                    particleEta[index]        = particle_list[h].eta;
                    particlePhi[index]        = particle_list[h].phi;
                    particlePhiDC[index]      = particle_list[h].phiDC;

                    particleErtTrigger[index] = particle_list[h].ertTrigger;

                    particleId[index]         = particle_list[h].id;

                    jetParticles_all.push_back(pseudoNeutral);
                    index++;
                }
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float chargedPt    = 0.0;
            float neutralPt    = 0.0;
            float discriminant = 0.0;

            float jetPt  = aFastJet.perp();
            float jetEta = aFastJet.pseudorapidity();
            float jetPhi = phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();

            bool ertTriggerFired = false;
            for (unsigned int iconst = 0; iconst < nconst; iconst++)
                {
                    unsigned int indx = constituents[iconst].user_index();
                    float deltaR      = dR(particleEta[indx], jetEta, particlePhi[indx], jetPhi);
                    float disc        = getDiscriminant(particlePt[indx], deltaR);
                    discriminant      += disc;

                    if (indx < indexCharged)//Charged particles
                        {
                            chargedPt += particlePt[indx];
                        }
                    else //Neutral particles
                        {
                            neutralPt += particlePt[indx];
                            if(particleErtTrigger[indx] == true)
                                {
                                    ertTriggerFired = true;
                                }
                        }
                }

            int jetArm = 1;
            if (jetPhi > 1.57)
                {
                    jetArm = 0;
                }

            //if(isData && !isMB && !ertTriggerFired)
	    //{
	    //continue;
	    //}

            float jetCf = chargedPt / jetPt;

            bool passJetLevelCuts = (jetPt > minPt) && ((float)nconst >= nc) && (jetCf >= minCf) && (jetCf <= maxCf);

            bool passAll = false;

            //Fiducial cuts: -1 is no fiducial cut, 0 is default, 1 is east arm, 2 is west arm, 3 is tight
            if(fiducialCut == -1)
                {
                    passAll = passJetLevelCuts;
                }
            if(fiducialCut == 0)
                {
                    passAll = passJetLevelCuts && passFiducialCut(jetEta, jetPhi);
                }
            if(fiducialCut == 1)
                {
                    passAll = passJetLevelCuts && passFiducialCut(jetEta, jetPhi) && (jetArm == 0);
                }
            if(fiducialCut == 2)
                {
                    passAll = passJetLevelCuts && passFiducialCut(jetEta, jetPhi) && (jetArm == 1);
                }
            if(fiducialCut == 3)
                {
                    passAll = passJetLevelCuts && passFiducialCutTight(jetEta, jetPhi);
                }

            if(passAll)
                {
                    jets temp;
                    temp.arm           = jetArm;
                    temp.centralityBin = -1;
                    temp.pT            = jetPt;
                    temp.eta           = jetEta;
                    temp.phi           = jetPhi;
                    temp.nc            = (float)nconst;
                    temp.cf            = jetCf;
                    temp.nf            = neutralPt / jetPt;
                    temp.disc          = discriminant;

                    antikt_jet_list.push_back(temp);

                    for (unsigned int iconst = 0; iconst < nconst; iconst++)
                        {
                            unsigned int indx = constituents[iconst].user_index();

                            particles tempConst;
                            tempConst.arm          = particleArm[indx];
                            tempConst.charge       = particleCharge[indx];

                            tempConst.energy       = particleEnergy[indx];
                            tempConst.mom          = particleMom[indx];
                            tempConst.pT           = particlePt[indx];
                            tempConst.eT           = particleEt[indx];
                            tempConst.px           = particlePx[indx];
                            tempConst.py           = particlePy[indx];
                            tempConst.pz           = particlePz[indx];
                            tempConst.eta          = particleEta[indx];
                            tempConst.phi          = particlePhi[indx];
                            tempConst.phiDC        = particlePhiDC[indx];

                            tempConst.ertTrigger   = particleErtTrigger[indx];

                            tempConst.id           = particleId[indx];
                            tempConst.jetPt        = jetPt;

                            constituent_list.push_back(tempConst);
                        }
                }
        }
    delete antikt;
}

void JetAnalyzer::GetAntiKt(std::vector<particles> particle_list,
                            fastjet::JetDefinition *antikt,
                            std::vector<jets>& antikt_jet_list,
                            TH1F *hAntikt[16],
                            TH1F *hAntiktCentrality[5],
                            TH1F *hAntiktErtEfficiency[5],
                            TH1F *hCharged[5],
                            TH1F *hNeutral[5])
{
    //Generate input for jet reconstuction
    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    //Sort the particle list- charged particles first and then neurtal
    std::sort(particle_list.begin(), particle_list.end(), sortParticle());

    unsigned int indexTotal = 0;
    unsigned int indexCharged = 0;
    if(particle_list.size() != 0)
        {
            for (unsigned int h = 0; h < particle_list.size(); h++)
                {
                    indexTotal++;
                    if (particle_list[h].charge != 0)
                        {
                            indexCharged++;
                        }
                }
        }
    else
        {
            indexTotal = 1;
        }

    float particlePt[indexTotal];
    fill(particlePt, particlePt + indexTotal / sizeof(float), -999.9);
    float particleEta[indexTotal];
    fill(particleEta, particleEta + indexTotal / sizeof(float), -999.9);
    float particlePhi[indexTotal];
    fill(particlePhi, particlePhi + indexTotal / sizeof(float), -999.9);
    bool particleErtTrigger[indexTotal];
    fill(particleErtTrigger, particleErtTrigger + indexTotal / sizeof(bool), false);

    int index = 0;
    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            if (particle_list[h].charge != 0)
                {
                    fastjet::PseudoJet pseudoCharged(particle_list[h].px,
                                                     particle_list[h].py,
                                                     particle_list[h].pz,
                                                     particle_list[h].mom);
                    pseudoCharged.set_user_index(index);
                    particlePt[index]           = particle_list[h].pT;
                    particleEta[index]          = particle_list[h].eta;
                    particlePhi[index]          = particle_list[h].phi;
                    particleErtTrigger[index]   = particle_list[h].ertTrigger;
                    jetParticles_all.push_back(pseudoCharged);
                    index++;
                }
            else
                {
                    fastjet::PseudoJet pseudoNeutral(particle_list[h].px,
                                                     particle_list[h].py,
                                                     particle_list[h].pz,
                                                     particle_list[h].energy);
                    pseudoNeutral.set_user_index(index);
                    particlePt[index]           = particle_list[h].pT;
                    particleEta[index]          = particle_list[h].eta;
                    particlePhi[index]          = particle_list[h].phi;
                    particleErtTrigger[index]   = particle_list[h].ertTrigger;
                    jetParticles_all.push_back(pseudoNeutral);
                    index++;
                }
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float chargedPt    = 0.0;
            float neutralPt    = 0.0;
            float discriminant = 0.0;

            float jetPt  = aFastJet.perp();
            float jetEta = aFastJet.pseudorapidity();
            float jetPhi = phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();

            bool ertTriggerFired = false;
            for (unsigned int iconst = 0; iconst < nconst; iconst++)
                {
                    unsigned int indx = constituents[iconst].user_index();
                    float deltaR      = dR(particleEta[indx], jetEta, particlePhi[indx], jetPhi);
                    float disc        = getDiscriminant(particlePt[indx], deltaR);
                    discriminant      += disc;

                    if (indx < indexCharged)//Charged particles
                        {
                            chargedPt += particlePt[indx];
                        }
                    else //Neutral particles
                        {
                            neutralPt += particlePt[indx];
                            if(particleErtTrigger[indx] == true)
                                {
                                    ertTriggerFired = true;
                                }
                        }
                }

            int arm = 1;
            if (jetPhi > 1.57)
                {
                    arm = 0;
                }

            jets temp;
            temp.arm           = arm;
            temp.centralityBin = centralityBin;
            temp.pT            = jetPt;
            temp.eta           = jetEta;
            temp.phi           = jetPhi;
            temp.nc            = (float)nconst;
            temp.cf            = chargedPt / jetPt;
            temp.nf            = neutralPt / jetPt;
            temp.disc          = discriminant;

            if(!isMB && !ertTriggerFired)
                {
                    continue;
                }

            bool jetLevelCut = temp.pT > 5.0 && temp.nc >= 3.0 && temp.cf >= 0.2 && temp.cf <= 0.7;
            bool passFiducial = passFiducialCut(jetEta, jetPhi);

            if(jetLevelCut && passFiducial)
                {
                    if(temp.pT > MINPT_RECO)
                        {
                            antikt_jet_list.push_back(temp);

                            hAntikt[1]->Fill(temp.disc);
                            hAntikt[2]->Fill(temp.eta);
                            hAntikt[3]->Fill(temp.phi);
                        }

                    //Centrality plots
                    hAntiktCentrality[0]->Fill(temp.pT);
                    if(ertTriggerFired)
                        {
                            hAntiktErtEfficiency[0]->Fill(temp.pT);
                        }
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hAntiktCentrality[centralityBin]->Fill(temp.pT);
                            if(ertTriggerFired)
                                {
                                    hAntiktErtEfficiency[centralityBin]->Fill(temp.pT);
                                }
                        }
                }

            if (temp.pT > 5.0 && temp.nc >= 3.0)
                {
                    hAntikt[0]->Fill(temp.disc);
                    hAntikt[4]->Fill(temp.cf);
                    hAntikt[6]->Fill(temp.nf);

                    if(temp.pT > MINPT_RECO)
                        {
                            hAntikt[5]->Fill(temp.cf);
                            hAntikt[7]->Fill(temp.nf);
                        }

                    hAntikt[8]->Fill(temp.pT);
                    if ((temp.cf > 0.1) && (temp.cf < 0.9))
                        {
                            hAntikt[9]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.8))
                        {
                            hAntikt[10]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.7))
                        {
                            hAntikt[11]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.6))
                        {
                            hAntikt[12]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.7) && (temp.disc > 15.0))
                        {
                            hAntikt[13]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.7) && (temp.disc > 20.0))
                        {
                            hAntikt[14]->Fill(temp.pT);
                        }
                    if ((temp.cf > 0.2) && (temp.cf < 0.7) && (temp.disc > 25.0))
                        {
                            hAntikt[15]->Fill(temp.pT);
                        }
                }

            if (temp.pT > MINPT_RECO && temp.nc >= 3.0)
                {
                    //Study the charged and neutral particles inside jet
                    for (unsigned int iconst = 0; iconst < nconst; iconst++)
                        {
                            unsigned int indx = constituents[iconst].user_index();
                            if (indx < indexCharged)//Charged particles
                                {
                                    hCharged[0]->Fill(particlePt[indx]);
                                    if ((temp.cf > 0.1) && (temp.cf < 0.9))
                                        {
                                            hCharged[1]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.8))
                                        {
                                            hCharged[2]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.7))
                                        {
                                            hCharged[3]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.6))
                                        {
                                            hCharged[4]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                }
                            else //Neutral particles
                                {
                                    hNeutral[0]->Fill(particlePt[indx]);
                                    if ((temp.cf > 0.1) && (temp.cf < 0.9))
                                        {
                                            hNeutral[1]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.8))
                                        {
                                            hNeutral[2]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.7))
                                        {
                                            hNeutral[3]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                    if ((temp.cf > 0.2) && (temp.cf < 0.6))
                                        {
                                            hNeutral[4]->Fill(particlePt[indx], 1 / particlePt[indx]);
                                        }
                                }
                        }
                }
        }
}

void JetAnalyzer::GetFakeJets(std::vector<particles> particle_list,
                              float R, float nc, float minPt, float minCf, float maxCf,
                              std::vector<particles>& shuffled_particle_list,
                              std::vector<jets>& fake_jet_list,
                              std::vector<particles>& constituent_list,
                              int fiducialCut)
{
    //Sort the particle list- charged particles first and then neurtal
    std::sort(particle_list.begin(), particle_list.end(), sortParticle());

    unsigned int indexChargedEast = 0;
    unsigned int indexChargedWest = 0;
    unsigned int indexNeutralEast = 0;
    unsigned int indexNeutralWest = 0;
    bool foundChargedEast = false;
    bool foundChargedWest = false;
    bool foundNeutralEast = false;
    bool foundNeutralWest = false;
    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            if (particle_list[h].charge != 0)
                {
                    if(particle_list[h].arm == 0)
                        {
                            indexChargedEast++;
                            foundChargedEast = true;
                        }
                    else
                        {
                            indexChargedWest++;
                            foundChargedWest = true;
                        }
                }
            else
                {
                    if(particle_list[h].arm == 0)
                        {
                            indexNeutralEast++;
                            foundNeutralEast = true;
                        }
                    else
                        {
                            indexNeutralWest++;
                            foundNeutralWest = true;
                        }
                }
        }

    if(!foundChargedEast)
        {
            indexChargedEast = 1;
        }
    if(!foundChargedWest)
        {
            indexChargedWest = 1;
        }

    if(!foundNeutralEast)
        {
            indexNeutralEast = 1;
        }
    if(!foundNeutralWest)
        {
            indexNeutralWest = 1;
        }

    float chargedParticleEtaEast[indexChargedEast];
    fill(chargedParticleEtaEast, chargedParticleEtaEast + indexChargedEast / sizeof(float), -999.9);
    float chargedParticlePhiEast[indexChargedEast];
    fill(chargedParticlePhiEast, chargedParticlePhiEast + indexChargedEast / sizeof(float), -999.9);

    float chargedParticleEtaWest[indexChargedWest];
    fill(chargedParticleEtaWest, chargedParticleEtaWest + indexChargedWest / sizeof(float), -999.9);
    float chargedParticlePhiWest[indexChargedWest];
    fill(chargedParticlePhiWest, chargedParticlePhiWest + indexChargedWest / sizeof(float), -999.9);

    float neutralParticleEtaEast[indexNeutralEast];
    fill(neutralParticleEtaEast, neutralParticleEtaEast + indexNeutralEast / sizeof(float), -999.9);
    float neutralParticlePhiEast[indexNeutralEast];
    fill(neutralParticlePhiEast, neutralParticlePhiEast + indexNeutralEast / sizeof(float), -999.9);

    float neutralParticleEtaWest[indexNeutralWest];
    fill(neutralParticleEtaWest, neutralParticleEtaWest + indexNeutralWest / sizeof(float), -999.9);
    float neutralParticlePhiWest[indexNeutralWest];
    fill(neutralParticlePhiWest, neutralParticlePhiWest + indexNeutralWest / sizeof(float), -999.9);

    unsigned int indexChargedNewEast = 0;
    unsigned int indexChargedNewWest = 0;
    unsigned int indexNeutralNewEast = 0;
    unsigned int indexNeutralNewWest = 0;
    for (unsigned int index = 0; index < particle_list.size(); index++)
        {
            if (particle_list[index].charge != 0)
                {
                    if (particle_list[index].arm == 0)
                        {
                            chargedParticleEtaEast[indexChargedNewEast] = particle_list[index].eta;
                            chargedParticlePhiEast[indexChargedNewEast] = particle_list[index].phi;
                            indexChargedNewEast++;
                        }
                    else
                        {
                            chargedParticleEtaWest[indexChargedNewWest] = particle_list[index].eta;
                            chargedParticlePhiWest[indexChargedNewWest] = particle_list[index].phi;
                            indexChargedNewWest++;
                        }
                }
            else
                {
                    if (particle_list[index].arm == 0)
                        {
                            neutralParticleEtaEast[indexNeutralNewEast] = particle_list[index].eta;
                            neutralParticlePhiEast[indexNeutralNewEast] = particle_list[index].phi;
                            indexNeutralNewEast++;
                        }
                    else
                        {
                            neutralParticleEtaWest[indexNeutralNewWest] = particle_list[index].eta;
                            neutralParticlePhiWest[indexNeutralNewWest] = particle_list[index].phi;
                            indexNeutralNewWest++;
                        }
                }
        }

    //Randomly swap (eta, phi) of tracks and clusters separately
    // Use a different seed value so that we don't get same result each time
    srand(time(NULL));

    //Start from the last element- don't need to run first element- hence i>0
    for (int i = indexChargedEast - 1; i > 0; i--)
        {
            // Pick a random index from 0 to i
            int indexRandomEast = rand() % (i + 1);

            //Swap (eta, phi) together
            float tempEtaEast = chargedParticleEtaEast[i];
            chargedParticleEtaEast[i] = chargedParticleEtaEast[indexRandomEast];
            chargedParticleEtaEast[indexRandomEast] = tempEtaEast;

            float tempPhiEast = chargedParticlePhiEast[i];
            chargedParticlePhiEast[i] = chargedParticlePhiEast[indexRandomEast];
            chargedParticlePhiEast[indexRandomEast] = tempPhiEast;
        }

    for (int i = indexChargedWest - 1; i > 0; i--)
        {
            int indexRandomWest = rand() % (i + 1);

            float tempEtaWest = chargedParticleEtaWest[i];
            chargedParticleEtaWest[i] = chargedParticleEtaWest[indexRandomWest];
            chargedParticleEtaWest[indexRandomWest] = tempEtaWest;

            float tempPhiWest = chargedParticlePhiWest[i];
            chargedParticlePhiWest[i] = chargedParticlePhiWest[indexRandomWest];
            chargedParticlePhiWest[indexRandomWest] = tempPhiWest;
        }

    for (int i = indexNeutralEast - 1; i > 0; i--)
        {
            int indexRandomEast = rand() % (i + 1);

            float tempEtaEast = neutralParticleEtaEast[i];
            neutralParticleEtaEast[i] = neutralParticleEtaEast[indexRandomEast];
            neutralParticleEtaEast[indexRandomEast] = tempEtaEast;

            float tempPhiEast = neutralParticlePhiEast[i];
            neutralParticlePhiEast[i] = neutralParticlePhiEast[indexRandomEast];
            neutralParticlePhiEast[indexRandomEast] = tempPhiEast;
        }

    for (int i = indexNeutralWest - 1; i > 0; i--)
        {
            int indexRandomWest = rand() % (i + 1);

            float tempEtaWest = neutralParticleEtaWest[i];
            neutralParticleEtaWest[i] = neutralParticleEtaWest[indexRandomWest];
            neutralParticleEtaWest[indexRandomWest] = tempEtaWest;

            float tempPhiWest = neutralParticlePhiWest[i];
            neutralParticlePhiWest[i] = neutralParticlePhiWest[indexRandomWest];
            neutralParticlePhiWest[indexRandomWest] = tempPhiWest;
        }

    //Anti-kt takes px, py, pz, mom/energy as input- so recalculate TLorentzVector with new eta and phi
    particles tempParticle;
    unsigned int iChargedEast = 0;
    unsigned int iChargedWest = 0;
    unsigned int iNeutralEast = 0;
    unsigned int iNeutralWest = 0;
    for(unsigned int iIndex = 0; iIndex < particle_list.size(); iIndex++)
        {
            tempParticle.arm        = particle_list[iIndex].arm;
            tempParticle.charge     = particle_list[iIndex].charge;
            tempParticle.phiDC      = particle_list[iIndex].phiDC;
            tempParticle.ertTrigger = particle_list[iIndex].ertTrigger;
            tempParticle.id         = particle_list[iIndex].id;
            tempParticle.jetPt      = particle_list[iIndex].jetPt;

            if(tempParticle.charge != 0)
                {
                    TLorentzVector *fVector = new TLorentzVector();
                    if(tempParticle.arm == 0)
                        {
                            fVector->SetPtEtaPhiE(particle_list[iIndex].mom * sin(2 * atan(exp(-chargedParticleEtaEast[iChargedEast]))),
                                                  chargedParticleEtaEast[iChargedEast],
                                                  chargedParticlePhiEast[iChargedEast],
                                                  particle_list[iIndex].mom);
                            iChargedEast++;
                        }
                    else
                        {
                            fVector->SetPtEtaPhiE(particle_list[iIndex].mom * sin(2 * atan(exp(-chargedParticleEtaWest[iChargedWest]))),
                                                  chargedParticleEtaWest[iChargedWest],
                                                  chargedParticlePhiWest[iChargedWest],
                                                  particle_list[iIndex].mom);
                            iChargedWest++;
                        }
                    tempParticle.energy     = particle_list[iIndex].energy;
                    tempParticle.mom        = fVector->P();
                    tempParticle.px         = fVector->Px();
                    tempParticle.py         = fVector->Py();
                    tempParticle.pz         = fVector->Pz();

                    tempParticle.phi        = phiReduce(fVector->Phi());
                    tempParticle.pT         = fVector->Pt();
                    tempParticle.eT         = particle_list[iIndex].energy * sin(fVector->Theta());
                    tempParticle.eta        = fVector->Eta();

                    delete fVector;
                }
            else
                {
                    TLorentzVector *fVector = new TLorentzVector();
                    if(tempParticle.arm == 0)
                        {
                            fVector->SetPtEtaPhiE(particle_list[iIndex].energy * sin(2 * atan(exp(-neutralParticleEtaEast[iNeutralEast]))),
                                                  neutralParticleEtaEast[iNeutralEast],
                                                  neutralParticlePhiEast[iNeutralEast],
                                                  particle_list[iIndex].energy);
                            iNeutralEast++;
                        }
                    else
                        {
                            fVector->SetPtEtaPhiE(particle_list[iIndex].energy * sin(2 * atan(exp(-neutralParticleEtaWest[iNeutralWest]))),
                                                  neutralParticleEtaWest[iNeutralWest],
                                                  neutralParticlePhiWest[iNeutralWest],
                                                  particle_list[iIndex].energy);
                            iNeutralWest++;
                        }
                    tempParticle.energy     = fVector->E();
                    tempParticle.mom        = fVector->P();
                    tempParticle.px         = fVector->Px();
                    tempParticle.py         = fVector->Py();
                    tempParticle.pz         = fVector->Pz();

                    tempParticle.phi        = phiReduce(fVector->Phi());
                    tempParticle.pT         = fVector->Pt();
                    tempParticle.eT         = fVector->Et();
                    tempParticle.eta        = fVector->Eta();

                    delete fVector;
                }
            shuffled_particle_list.push_back(tempParticle);
        }

    //Reconstruct Anti-kt jet
    GetAntiKtCommon(shuffled_particle_list, R, nc, minPt, minCf, maxCf, fake_jet_list, constituent_list,
                    fiducialCut);
}

void JetAnalyzer::GetPi0(std::vector<clusters> cluster_list,
                         TH1F *hMass[5][8],
                         TH2F *hMassVsPt[5][8],
                         TH1F *hAsymmetry[5][8])
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
                    tempPhoton.ertTrigger = cluster_list[c].ert4x4Bit;

                    all_photons.push_back(tempPhoton);
                }
        }

    //Make pi0s
    for (unsigned int p = 0; p < all_photons.size(); p++)
        {
            for (unsigned int pp = p + 1; pp < all_photons.size(); pp++)
                {
                    //Require same armsect
                    if(all_photons[p].armsect != all_photons[pp].armsect)
                        {
                            continue;
                        }

                    //For ERT dataset, require atleast one photon to fire ERT
                    bool eFired = all_photons[p].ertTrigger || all_photons[pp].ertTrigger;
                    if(!isMB && !eFired)
                        {
                            continue;
                        }

                    float asymmetry = fabs((all_photons[p].pT - all_photons[pp].pT) / (all_photons[p].pT + all_photons[pp].pT));
                    hAsymmetry[0][all_photons[p].armsect]->Fill(asymmetry);
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hAsymmetry[centralityBin][all_photons[p].armsect]->Fill(asymmetry);
                        }

                    if(asymmetry > 0.8)
                        {
                            continue;
                        }

                    TLorentzVector photon1;
                    photon1.SetPtEtaPhiM(all_photons[p].pT, all_photons[p].eta, all_photons[p].phi, 0.0);
                    TLorentzVector photon2;
                    photon2.SetPtEtaPhiM(all_photons[pp].pT, all_photons[pp].eta, all_photons[pp].phi, 0.0);

                    TLorentzVector pion = photon1 + photon2;

                    float pionMass = pion.M();
                    float pionPt = pion.Pt();
                    hMassVsPt[0][all_photons[p].armsect]->Fill(pionMass, pionPt);
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hMassVsPt[centralityBin][all_photons[p].armsect]->Fill(pionMass, pionPt);
                        }

                    if(pionPt < 3.0 || pionPt > 10.0)
                        {
                            continue;
                        }

                    hMass[0][all_photons[p].armsect]->Fill(pionMass);
                    if(centralityBin >= 1 && centralityBin <= 4)
                        {
                            hMass[centralityBin][all_photons[p].armsect]->Fill(pionMass);
                        }
                }
        }
}


void JetAnalyzer::InitTrees(bool writeTrees)
{
    if(writeTrees)
        {
            //TTree for tracks and clusters
            allTracks = new TTree("allTracks", "Central arm tracks");
            allTracks->Branch("allTracks", &all_tracks);

            allCharged = new TTree("allCharged", "Charged particles");
            allCharged->Branch("allCharged", &charged_particles);

            allClusters = new TTree("allClusters", "Central arm clusters");
            allClusters->Branch("allClusters", &all_clusters);

            allNeutrals = new TTree("allNeutrals", "Neutral particles");
            allNeutrals->Branch("allCharged", &neutral_particles);

            allparticles = new TTree("allparticles", "Neutral+Charged particles");
            allparticles->Branch("allparticles", &all_particles);

            //TTree for aniti-kt jets
            antikt0 = new TTree("antikt0", "Jets, anti-kt, R = 0.15");
            antikt0->Branch("antikt0", &antikt_jets_0);

            antikt1 = new TTree("antikt1", "Jets, anti-kt, R = 0.2");
            antikt1->Branch("antikt1", &antikt_jets_1);

            antikt2 = new TTree("antikt2", "Jets, anti-kt, R = 0.25");
            antikt2->Branch("antikt2", &antikt_jets_2);

            antikt3 = new TTree("antikt3", "Jets, anti-kt, R = 0.3");
            antikt3->Branch("antikt3", &antikt_jets_3);
        }
}

void JetAnalyzer::FillTrees(bool writeTrees)
{
    if(writeTrees)
        {
            allTracks->Fill();
            allClusters->Fill();
            allCharged->Fill();
            allNeutrals->Fill();

            allparticles->Fill();

            antikt0->Fill();
            antikt1->Fill();
            antikt2->Fill();
            antikt3->Fill();
        }
}

void JetAnalyzer::InitHistograms()
{
    //Histograms
    //General Histos
    hEvents = new TH1F("hEvents", "Number of events", 100, 0, 100);

    hEffectiveEvents = new TH1F("hEffectiveEvents", "Effective number of events (p+p, MB)", 5, 0, 5);

    //For Pair cut
    hAllPairs = new TH2F("hAllPairs", "All pairs", 30, -0.5, 0.5, 100, -0.1, 0.1);
    hSameChargePairs = new TH2F("hSameChargePairs", "Same Charge pairs", 30, -0.5, 0.5, 100, -0.1, 0.1);
    hOppositeChargePairs = new TH2F("hOppositeChargePairs", "Opposite Charge pairs", 30, -0.5, 0.5, 100, -0.1, 0.1);

    //For Matching Study
    hEmcsdPhi = new TH2F("hEmcsdPhi", "EMCal sdPhi", 20, -5.00, 5.00, 500, 0.00, 25.00);
    hEmcsdZ = new TH2F("hEmcsdZ", "EMCal sdZ", 20, -5.00, 5.00, 500, 0.00, 25.00);

    hPc3sdPhi = new TH2F("hPc3sdPhi", "EMCal sdPhi", 20, -5.00, 5.00, 500, 0.00, 25.00);
    hPc3sdZ = new TH2F("hPc3sdZ", "EMCal sdZ", 20, -5.00, 5.00, 500, 0.00, 25.00);

    //To make quality mask
    hMakeMask_NE_X1 = new TH2F("hMakeMask_NE_X1", "Tracks with X1 bit, NE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SE_X1 = new TH2F("hMakeMask_SE_X1", "Tracks with X1 bit, SE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_NW_X1 = new TH2F("hMakeMask_NW_X1", "Tracks with X1 bit, NW",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SW_X1 = new TH2F("hMakeMask_SW_X1", "Tracks with X1 bit, SW",  400, 0, 80, 120, -0.6, 0.6);

    hMakeMask_NE_X2 = new TH2F("hMakeMask_NE_X2", "Tracks with X2 bit, NE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SE_X2 = new TH2F("hMakeMask_SE_X2", "Tracks with X2 bit, SE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_NW_X2 = new TH2F("hMakeMask_NW_X2", "Tracks with X2 bit, NW",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SW_X2 = new TH2F("hMakeMask_SW_X2", "Tracks with X2 bit, SW",  400, 0, 80, 120, -0.6, 0.6);

    hMakeMask_NE_UV = new TH2F("hMakeMask_NE_UV", "Tracks with UV bit, NE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SE_UV = new TH2F("hMakeMask_SE_UV", "Tracks with UV bit, SE",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_NW_UV = new TH2F("hMakeMask_NW_UV", "Tracks with UV bit, NW",  400, 0, 80, 120, -0.6, 0.6);
    hMakeMask_SW_UV = new TH2F("hMakeMask_SW_UV", "Tracks with UV bit, SW",  400, 0, 80, 120, -0.6, 0.6);

    //To check quality mask
    hCheckMask_NE_X1 = new TH2F("hCheckMask_NE_X1", "Tracks not in broken X1, NE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SE_X1 = new TH2F("hCheckMask_SE_X1", "Tracks not in broken X1, SE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_NW_X1 = new TH2F("hCheckMask_NW_X1", "Tracks not in broken X1, NW",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SW_X1 = new TH2F("hCheckMask_SW_X1", "Tracks not in broken X1, SW",  400, 0, 80, 120, -0.6, 0.6);

    hCheckMask_NE_X2 = new TH2F("hCheckMask_NE_X2", "Tracks not in broken X2, NE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SE_X2 = new TH2F("hCheckMask_SE_X2", "Tracks not in broken X2, SE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_NW_X2 = new TH2F("hCheckMask_NW_X2", "Tracks not in broken X2, NW",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SW_X2 = new TH2F("hCheckMask_SW_X2", "Tracks not in broken X2, SW",  400, 0, 80, 120, -0.6, 0.6);

    hCheckMask_NE_UV = new TH2F("hCheckMask_NE_UV", "Tracks not in broken UV, NE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SE_UV = new TH2F("hCheckMask_SE_UV", "Tracks not in broken UV, SE",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_NW_UV = new TH2F("hCheckMask_NW_UV", "Tracks not in broken UV, NW",  400, 0, 80, 120, -0.6, 0.6);
    hCheckMask_SW_UV = new TH2F("hCheckMask_SW_UV", "Tracks not in broken UV, SW",  400, 0, 80, 120, -0.6, 0.6);

    //To compare
    hRegularQuality_NE = new TH2F("hRegularQuality_NE", "Quality 63||31, NE",  400, 0, 80, 120, -0.6, 0.6);
    hRegularQuality_SE = new TH2F("hRegularQuality_SE", "Quality 63||31, SE",  400, 0, 80, 120, -0.6, 0.6);
    hRegularQuality_NW = new TH2F("hRegularQuality_NW", "Quality 63||31, NW",  400, 0, 80, 120, -0.6, 0.6);
    hRegularQuality_SW = new TH2F("hRegularQuality_SW", "Quality 63||31, SW",  400, 0, 80, 120, -0.6, 0.6);

    hModifiedQuality_NE = new TH2F("hModifiedQuality_NE", "Modified Quality Cut, NE",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_SE = new TH2F("hModifiedQuality_SE", "Modified Quality Cut, SE",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_NW = new TH2F("hModifiedQuality_NW", "Modified Quality Cut, NW",  400, 0, 80, 120, -0.6, 0.6);
    hModifiedQuality_SW = new TH2F("hModifiedQuality_SW", "Modified Quality Cut, SW",  400, 0, 80, 120, -0.6, 0.6);

    //Tracks and clusters
    hCharged_ecore = new TH1F("hCharged_ecore", "Ecore of good charged tracks", 100, 0, 2);
    hCharged_EP1 = new TH1F("hCharged_EP1", "E/P for all good tracks", 100, 0, 2);
    hCharged_EP2 = new TH1F("hCharged_EP2", "E/P for pT<4.5 and n0>=2", 100, 0, 2);

    for (int ias = 0; ias < NARMSECT; ias++)
        {
            int ny, nz;
            if (EmcMap::IsPbGl(ias))
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
    hTowerE = new TH2F("hTowerE", "Ecore vs Tower-id", NTOWER, 0, NTOWER, 10, 0, 10.0);


    for (unsigned int c = 0; c <= 4; c++)
        {
            hTracksZed[c] = new TH1F(Form("hTracksZed_%u", c),
                                     "Tracks going into Jet reconstruction", 400, -100.0, 100.0);
            hTracksPhi[c] = new TH1F(Form("hTracksPhi_%u", c),
                                     "Tracks going into Jet reconstruction", 200, -1.0, 4.0);

            hClustersEta[c] = new TH1F(Form("hClustersEta_%u", c),
                                       "Clusters going into Jet reconstruction", 160, -0.4, 0.4);
            hClustersPhi[c] = new TH1F(Form("hClustersPhi_%u", c),
                                       "Clusters going into Jet reconstruction", 200, -1.0, 4.0);
            hClustersTower[c] = new TH1F(Form("hClustersTower_%u", c),
                                         "Hits vs Tower-id", NTOWER, 0, NTOWER);
        }


    //Jets
    for (unsigned int c = 0; c < 16; c++)
        {
            int bins = 0;
            float binMin = 0.0;
            float binMax = 0.0;
            char Title[256];

            if (c == 0)
                {
                    bins = 1000;
                    binMin = 0.0;
                    binMax = 100.0;
                    sprintf(Title, "Discriminant of Jet (nc>=3 && pT>5)");
                }
            if (c == 1)
                {
                    bins = 1000;
                    binMin = 0.0;
                    binMax = 100.0;
                    sprintf(Title, "Discriminant of Jet (Jet Level)");
                }
            if (c == 2)
                {
                    bins = 160;
                    binMin = -0.4;
                    binMax = 0.4;
                    sprintf(Title, "Eta of Jet (Jet Level)");
                }
            if (c == 3)
                {
                    bins = 200;
                    binMin = -1.0;
                    binMax = 4.0;
                    sprintf(Title, "Phi of Jet (Jet Level)");
                }
            if (c == 4)
                {
                    bins = 100;
                    binMin = -0.02;
                    binMax = 1.02;
                    sprintf(Title, "Charged Fraction of Jet (nc>=3 && pT>5)");
                }
            if (c == 5)
                {
                    bins = 100;
                    binMin = -0.02;
                    binMax = 1.02;
                    sprintf(Title, "Charged Fraction of Jet (nc>=3 && pT>15)");
                }
            if (c == 6)
                {
                    bins = 100;
                    binMin = -0.02;
                    binMax = 1.02;
                    sprintf(Title, "Neutral Fraction of Jet (nc>=3 && pT>5)");
                }
            if (c == 7)
                {
                    bins = 100;
                    binMin = -0.02;
                    binMax = 1.02;
                    sprintf(Title, "Neutral Fraction of Jet (nc>=3 && pT>15)");
                }
            if (c > 7)
                {
                    sprintf(Title, "p_{T} of Jet (cut %u)", c - 8);

                    hAntikt0[c] = new TH1F(Form("hAntikt_%d_%u", 15, c), Title, NPTBINS, PTBINS);
                    hAntikt1[c] = new TH1F(Form("hAntikt_%d_%u", 2, c), Title, NPTBINS, PTBINS);
                    hAntikt2[c] = new TH1F(Form("hAntikt_%d_%u", 25, c), Title, NPTBINS, PTBINS);
                    hAntikt3[c] = new TH1F(Form("hAntikt_%d_%u", 3, c), Title, NPTBINS, PTBINS);
                }

            if(c <= 7)
                {
                    hAntikt0[c] = new TH1F(Form("hAntikt_%d_%u", 15, c), Title, bins, binMin, binMax);
                    hAntikt1[c] = new TH1F(Form("hAntikt_%d_%u", 2, c), Title, bins, binMin, binMax);
                    hAntikt2[c] = new TH1F(Form("hAntikt_%d_%u", 25, c), Title, bins, binMin, binMax);
                    hAntikt3[c] = new TH1F(Form("hAntikt_%d_%u", 3, c), Title, bins, binMin, binMax);
                }
        }

    //Centrality
    for (unsigned int c = 0; c < 5; c++)
        {
            char TitleCent[256];
            sprintf(TitleCent, "p_{T} of Jet (Centrality: %u)", c);

            hAntiktCentrality0[c] = new TH1F(Form("hAntiktCentrality_%d_%u", 15, c), TitleCent, NPTBINS, PTBINS);
            hAntiktCentrality1[c] = new TH1F(Form("hAntiktCentrality_%d_%u", 2, c), TitleCent, NPTBINS, PTBINS);
            hAntiktCentrality2[c] = new TH1F(Form("hAntiktCentrality_%d_%u", 25, c), TitleCent, NPTBINS, PTBINS);
            hAntiktCentrality3[c] = new TH1F(Form("hAntiktCentrality_%d_%u", 3, c), TitleCent, NPTBINS, PTBINS);

            //ERT Efficiency
            char TitleEff[256];
            sprintf(TitleEff, "p_{T} of Jet (For Efficiency: %u)", c);

            hAntiktErtEfficiency0[c] = new TH1F(Form("hAntiktErtEfficiency_%d_%u", 15, c), TitleEff, NPTBINS, PTBINS);
            hAntiktErtEfficiency1[c] = new TH1F(Form("hAntiktErtEfficiency_%d_%u", 2, c), TitleEff, NPTBINS, PTBINS);
            hAntiktErtEfficiency2[c] = new TH1F(Form("hAntiktErtEfficiency_%d_%u", 25, c), TitleEff, NPTBINS, PTBINS);
            hAntiktErtEfficiency3[c] = new TH1F(Form("hAntiktErtEfficiency_%d_%u", 3, c), TitleEff, NPTBINS, PTBINS);
        }

    //Constituents study
    for (unsigned int c = 0; c < 5; c++)
        {
            hCharged0[c] = new TH1F(Form("hCharged_%d_%u", 15, c),
                                    Form("p_{T} of Charged Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hCharged1[c] = new TH1F(Form("hCharged_%d_%u", 2, c),
                                    Form("p_{T} of Charged Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hCharged2[c] = new TH1F(Form("hCharged_%d_%u", 25, c),
                                    Form("p_{T} of Charged Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hCharged3[c] = new TH1F(Form("hCharged_%d_%u", 3, c),
                                    Form("p_{T} of Charged Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);

            hNeutral0[c] = new TH1F(Form("hNeutral_%d_%u", 15, c),
                                    Form("p_{T} of Neutral Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hNeutral1[c] = new TH1F(Form("hNeutral_%d_%u", 2, c),
                                    Form("p_{T} of Neutral Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hNeutral2[c] = new TH1F(Form("hNeutral_%d_%u", 25, c),
                                    Form("p_{T} of Neutral Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
            hNeutral3[c] = new TH1F(Form("hNeutral_%d_%u", 3, c),
                                    Form("p_{T} of Neutral Particles inside Anti-kt Jet (cut %u)", c),
                                    25 * 4, 0, 25);
        }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Final Jets (Anti-kt, R=0.2)
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (unsigned int c = 0; c < 5; c++)
        {
            hFinalJets[c] = new TH1F(Form("hFinalJets_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsEtaPhi[c] = new TH2F(Form("hFinalJetsEtaPhi_%u", c), "", 160, -0.4, 0.4, 200, -1.0, 4.0);

            hFinalJetsDefault[c] = new TH1F(Form("hFinalJetsDefault_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsEast[c] = new TH1F(Form("hFinalJetsEast_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsWest[c] = new TH1F(Form("hFinalJetsWest_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsFidTight[c] = new TH1F(Form("hFinalJetsFidTight_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsNc[c] = new TH1F(Form("hFinalJetsNc_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsCf[c] = new TH1F(Form("hFinalJetsCf_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsNcCf[c] = new TH1F(Form("hFinalJetsNcCf_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFinalJetsTrClTight[c] = new TH1F(Form("hFinalJetsTrClTight_%u", c), "", NPTBINS_RECO, PTBINS_RECO);

            hFakeJetsDefault[c] = new TH1F(Form("hFakeJetsDefault_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsEast[c] = new TH1F(Form("hFakeJetsEast_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsWest[c] = new TH1F(Form("hFakeJetsWest_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsFidTight[c] = new TH1F(Form("hFakeJetsFidTight_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsNc[c] = new TH1F(Form("hFakeJetsNc_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsCf[c] = new TH1F(Form("hFakeJetsCf_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsNcCf[c] = new TH1F(Form("hFakeJetsNcCf_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
            hFakeJetsTrClTight[c] = new TH1F(Form("hFakeJetsTrClTight_%u", c), "", NPTBINS_RECO, PTBINS_RECO);
        }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Post QM15
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Vertex
    for (unsigned int c = 0; c < 5; c++)
        {
            hVertex[c] = new TH1F(Form("hVertex_%u", c), "Vertex Distribution, centrality", 100 * 2, -50, 50);
            hTotalVertexVsRun[c] = new TH1F(Form("hTotalVertexVsRun_%u", c), "Total Vertex vs Run, centrality", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
            hVertexEventsVsRun[c] = new TH1F(Form("hVertexEventsVsRun_%u", c), "Vertex events vs Run, centrality", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
        }

    //Centrality
    hCentrality = new TH1F("hCentrality", "Centrality Distribution", 120, -10, 110);
    hTotalCentralityVsRun = new TH1F("hTotalCentralityVsRun", "Total Centrality vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
    hCentralityEventsVsRun = new TH1F("hCentralityEventsVsRun", "Centrality events vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);

    hCentralityAfterVertex = new TH1F("hCentralityAfterVertex", "Centrality Distribution after vertex cut", 120, -10, 110);
    hTotalCentralityAfterVertexVsRun = new TH1F("hTotalCentralityAfterVertexVsRun", "Total Centrality Distribution after vertex cut vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);

    hEffectiveEventsVsRun = new TH1F("hEffectiveEventsVsRun", "Effective events vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);

    //Jet yield
    for (unsigned int c = 0; c < 5; c++)
        {
            hTotalEventsVsRun[c] = new TH1F(Form("hTotalEventsVsRun_%u", c), "Total Events vs Run, centrality", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
            hGoodEventsVsRun[c] = new TH1F(Form("hGoodEventsVsRun_%u", c), "Good Events vs Run, centrality", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);

            hJetYieldVsRun[c] = new TH1F(Form("hJetYieldVsRun_%u", c), "Jet Yield vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
            hJetYieldHighPtVsRun[c] = new TH1F(Form("hJetYieldHighPtVsRun_%u", c), "Jet Yield (pT>13.8) vs Run", RUN_BIN, IRUN_NUMBER, FRUN_NUMBER);
        }

    //Energy scale
    for (unsigned int c = 0; c < 5; c++)
        {
            for (unsigned int s = 0; s < 8; s++)
                {
                    //1. pi0 in EMCal
                    hPi0Mass[c][s] = new TH1F(Form("hPi0Mass_%u_%u", c, s), "Pi0 mass", 720, 0.05, 0.35);
                    hPi0MassVsPt[c][s] = new TH2F(Form("hPi0MassVsPt_%u_%u", c, s), "Pi0 mass vs Pt", 720, 0.05, 0.35, 100, 0, 25);
                    hPi0Asymmetry[c][s] = new TH1F(Form("hPi0Asymmetry_%u_%u", c, s), "Pt assymetry", 120, -0.1, 1.1);
                }
        }

    //R=0.3
    for (unsigned int c = 0; c < 5; c++)
        {
            hR3Jets[c] = new TH1F(Form("hR3Jets_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
            hR3FakeJets[c] = new TH1F(Form("hR3FakeJets_%u", c), "", NPTBINS_RECO_R3, PTBINS_RECO_R3);

            hR3AllJets[c] = new TH1F(Form("hR3AllJets_%u", c), "", NPTBINS, PTBINS);
            hR3ERTJets[c] = new TH1F(Form("hR3ERTJets_%u", c), "", NPTBINS, PTBINS);
        }
    hR3JetsEast = new TH1F("hR3JetsEast", "", NPTBINS_RECO_R3, PTBINS_RECO_R3);
    hR3JetsWest = new TH1F("hR3JetsWest", "", NPTBINS_RECO_R3, PTBINS_RECO_R3);

    //UE study
    for (unsigned int c = 0; c < 5; c++)
        {
            hBackground[c] = new TH1F(Form("hBackground_%u", c), "", 120, 0.0, 30.0);
        }
}

void JetAnalyzer::FillHistograms()
{
    for (unsigned int t = 0; t < all_tracks.size(); t++)
        {
            int quality         = all_tracks[t].quality;
            int arm             = all_tracks[t].arm;

            float zedDC         = all_tracks[t].zedDC;
            float alpha         = all_tracks[t].alpha;
            float pT            = all_tracks[t].pT;
            float board         = all_tracks[t].board;

            float emcsdz        = all_tracks[t].emcsdz;
            float emcsdphi      = all_tracks[t].emcsdphi;

            float pc3sdz        = all_tracks[t].pc3sdz;
            float pc3sdphi      = all_tracks[t].pc3sdphi;

            bool inBrokenX1     = all_tracks[t].inBrokenX1;
            bool inBrokenX2     = all_tracks[t].inBrokenX2;
            bool inBrokenUV     = all_tracks[t].inBrokenUV;
            bool passQuality    = all_tracks[t].passQuality;
            bool x1Used         = all_tracks[t].x1Used;
            bool x2Used         = all_tracks[t].x2Used;
            bool uvUsed         = all_tracks[t].uvUsed;

            bool passDC         = all_tracks[t].passDC;

            if(passDC)
                {
                    //For Matching study
                    hEmcsdPhi->Fill(emcsdphi, pT);
                    hEmcsdZ->Fill(emcsdz, pT);

                    hPc3sdPhi->Fill(pc3sdphi, pT);
                    hPc3sdZ->Fill(pc3sdz, pT);
                }


            //For Quality Mask
            if (pT > 0.2 && pT < 25)
                {
                    if (x1Used)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hMakeMask_NE_X1->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hMakeMask_NW_X1->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hMakeMask_SE_X1->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hMakeMask_SW_X1->Fill(board, alpha);
                                }
                        }

                    if (x2Used)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hMakeMask_NE_X2->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hMakeMask_NW_X2->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hMakeMask_SE_X2->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hMakeMask_SW_X2->Fill(board, alpha);
                                }
                        }
                    if (uvUsed)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hMakeMask_NE_UV->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hMakeMask_NW_UV->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hMakeMask_SE_UV->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hMakeMask_SW_UV->Fill(board, alpha);
                                }
                        }


                    if (!inBrokenX1)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hCheckMask_NE_X1->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hCheckMask_NW_X1->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hCheckMask_SE_X1->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hCheckMask_SW_X1->Fill(board, alpha);
                                }
                        }

                    if (!inBrokenX2)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hCheckMask_NE_X2->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hCheckMask_NW_X2->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hCheckMask_SE_X2->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hCheckMask_SW_X2->Fill(board, alpha);
                                }
                        }
                    if (!inBrokenUV)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hCheckMask_NE_UV->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hCheckMask_NW_UV->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hCheckMask_SE_UV->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hCheckMask_SW_UV->Fill(board, alpha);
                                }
                        }

                    if (quality == 63 || quality == 31)
                        {
                            if (zedDC > 0 && arm == 0)
                                {
                                    hRegularQuality_NE->Fill(board, alpha);
                                }
                            if (zedDC > 0 && arm == 1)
                                {
                                    hRegularQuality_NW->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 0)
                                {
                                    hRegularQuality_SE->Fill(board, alpha);
                                }
                            if (zedDC < 0 && arm == 1)
                                {
                                    hRegularQuality_SW->Fill(board, alpha);
                                }
                        }

                    if (passQuality)
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
        }

    //Tracks
    for (unsigned int t = 0; t < all_tracks.size(); t++)
        {
            int n0              = all_tracks[t].n0;
            float mom           = all_tracks[t].mom;
            float energy        = all_tracks[t].energy;
            float pT            = all_tracks[t].pT;
            float emcsdz        = all_tracks[t].emcsdz;
            float emcsdphi      = all_tracks[t].emcsdphi;
            bool passFirst      = all_tracks[t].passFirst;

            if(passFirst)
                {
                    bool emcMatching = sqrt((emcsdphi * emcsdphi) + (emcsdz * emcsdz)) < 3.0;

                    if (emcMatching && n0 <= 0)
                        {
                            hCharged_ecore->Fill(energy);
                        }
                    hCharged_EP1->Fill(energy / mom);
                    if (pT < 4.5 && n0 >= 2)
                        {
                            hCharged_EP2->Fill(energy / mom);
                        }
                }


            //For pair cut
            if (!all_tracks[t].passFirst)
                {
                    continue;
                }
            for (unsigned int tt = 0; tt < all_tracks.size(); tt++)
                {
                    if (!all_tracks[tt].passFirst)
                        {
                            continue;
                        }
                    if ((t == tt) || (all_tracks[t].arm != all_tracks[tt].arm))
                        {
                            continue;
                        }

                    float dPhi = all_tracks[t].phiDC - all_tracks[tt].phiDC;
                    float dZed = all_tracks[t].zedDC - all_tracks[tt].zedDC;

                    hAllPairs->Fill(dZed, dPhi);

                    if(all_tracks[t].charge == all_tracks[tt].charge)
                        {
                            hSameChargePairs->Fill(dZed, dPhi);
                        }

                    if(all_tracks[t].charge != all_tracks[tt].charge)
                        {
                            hOppositeChargePairs->Fill(dZed, dPhi);
                        }
                }
        }

    for (unsigned int t = 0; t < charged_particles.size(); t++)
        {
            float zedDC           = charged_particles[t].zedDC;
            float phiDC           = charged_particles[t].phiDC;

            hTracksZed[0]->Fill(zedDC);
            hTracksPhi[0]->Fill(phiDC);
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    hTracksZed[centralityBin]->Fill(zedDC);
                    hTracksPhi[centralityBin]->Fill(phiDC);
                }
        }

    //Clusters
    for (int iclus = 0; iclus < all_clusters.size(); iclus++)
        {
            if(all_clusters[iclus].passEverything)
                {
                    int armsect      = all_clusters[iclus].armsect;
                    int yTowerPos    = all_clusters[iclus].yTowerPos;
                    int zTowerPos    = all_clusters[iclus].zTowerPos;
                    hSectorHits[armsect]->Fill(zTowerPos, yTowerPos);

                    float energy = all_clusters[iclus].energy;
                    int towerId = all_clusters[iclus].towerId;

                    hTowerE->Fill(towerId, energy);
                }
        }

    for (unsigned int c = 0; c < neutral_particles.size(); c++)
        {
            float eta = neutral_particles[c].eta;
            float phi = neutral_particles[c].phi;
            int towerId = neutral_particles[c].towerId;

            hClustersEta[0]->Fill(eta);
            hClustersPhi[0]->Fill(phi);
            hClustersTower[0]->Fill(towerId);
            if(centralityBin >= 1 && centralityBin <= 4)
                {
                    hClustersEta[centralityBin]->Fill(eta);
                    hClustersPhi[centralityBin]->Fill(phi);
                    hClustersTower[centralityBin]->Fill(towerId);
                }
        }
}


















































