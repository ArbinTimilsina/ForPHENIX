//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//HepMc tools
#include <PHHepMCGenEvent.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//PHPythia tools
#include <PHPythiaHeader.h>
#include <PHPythiaContainerV2.h>
#include <PHPyCommon.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//C tools
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <algorithm> // remove and remove_if

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TNtuple.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//Gaussina Filter tool
#include <jetevent/event.h>

//My source file
#include "JetSimWithoutDetector.h"

using namespace std;
using namespace findNode;
using namespace fastjet;

//================================ Constructor ================================
//Here we can initiate some variables
JetSimWithoutDetector::JetSimWithoutDetector(string outfilename)
    : SubsysReco("JetSimWithoutDetector"),
      verbo(1),
      outfname(outfilename),
      fillTrees(false)
{
    return;
}

int JetSimWithoutDetector::Init(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    outfile = new TFile(outfname.c_str(), "RECREATE");

    //For information on primary
    accept = new TAcceptParticle();
    tpythia6 = new TPythia6();

    //Jet algorithm definition
    //Anit-kt
    antikt_00 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.15, fastjet::E_scheme, fastjet::Best);
    antikt_01 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.2, fastjet::E_scheme, fastjet::Best);
    antikt_02 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.25, fastjet::E_scheme, fastjet::Best);
    antikt_03 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.3, fastjet::E_scheme, fastjet::Best);

    //Gaussian Filter
    background_model_perp = new jet::factorized_background_model_t();
    background_model_time = new jet::factorized_background_model_t();
    background_model_z = new jet::factorized_background_model_t();

    static const std::pair<float, float> phenix_acceptance_pseudorapidity_range = std::pair<float, float>(-8.0F, 8.0F);
    static const int phenix_acceptance_npixel_pseudorapidity = 640;

    filter_00 = new jet::reconstruction_filtering_iir_t(0.106F, phenix_acceptance_pseudorapidity_range,
							phenix_acceptance_npixel_pseudorapidity,
							*background_model_perp,
							*background_model_time,
							*background_model_z);

    filter_01 = new jet::reconstruction_filtering_iir_t(0.141F, phenix_acceptance_pseudorapidity_range,
							phenix_acceptance_npixel_pseudorapidity,
							*background_model_perp,
							*background_model_time,
							*background_model_z);

    filter_02 = new jet::reconstruction_filtering_iir_t(0.177F, phenix_acceptance_pseudorapidity_range,
							phenix_acceptance_npixel_pseudorapidity,
							*background_model_perp,
							*background_model_time,
							*background_model_z);

    filter_03 = new jet::reconstruction_filtering_iir_t(0.212F, phenix_acceptance_pseudorapidity_range,
							phenix_acceptance_npixel_pseudorapidity,
							*background_model_perp,
							*background_model_time,
							*background_model_z);

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    //General
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    hForCentrality = new TH1D("hForCentrality", "sHIJING Final Particle Multiplicity", 200, 0, 8000);
    hImpactParameter = new TH1D("hImpactParameter", "sHIJING Impact Parameter", 210, -1, 20);
    hEventPlane = new TH1D("hEventPlane", "sHIJING Event Plane", 140, -0.5, 6.5);

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //sHijing
    //*****************************************************************************************************************************
    hDistance_hijing_Antikt_0 = new TH1D("hDistance_hijing_Antikt_0",
                                         "#DeltaR between sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True sHIJING Jet. Anti-kt, R=0.15",
                                         600, 0.0, 0.6);
    hDistance_hijing_Antikt_1 = new TH1D("hDistance_hijing_Antikt_1",
                                         "#DeltaR between sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True sHIJING Jet. Anti-kt, R=0.2",
                                         600, 0.0, 0.6);
    hDistance_hijing_Antikt_2 = new TH1D("hDistance_hijing_Antikt_2",
                                         "#DeltaR between sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True sHIJING Jet. Anti-kt, R=0.25",
                                         600, 0.0, 0.6);
    hDistance_hijing_Antikt_3 = new TH1D("hDistance_hijing_Antikt_3",
                                         "#DeltaR between sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True sHIJING Jet. Anti-kt, R=0.3",
                                         600, 0.0, 0.6);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtTrueJet_hijing_Antikt_0 = new TH1D("hPtTrueJet_hijing_Antikt_0",
                                          "p_{T} of sHIJING True Jet, Antikt, R=0.15; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_hijing_Antikt_1 = new TH1D("hPtTrueJet_hijing_Antikt_1",
                                          "p_{T} of sHIJING True Jet, Antikt, R=0.2; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_hijing_Antikt_2 = new TH1D("hPtTrueJet_hijing_Antikt_2",
                                          "p_{T} of sHIJING True Jet, Antikt, R=0.25; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_hijing_Antikt_3 = new TH1D("hPtTrueJet_hijing_Antikt_3",
                                          "p_{T} of sHIJING True Jet, Antikt, R=0.3; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtMatchedJet_hijing_Antikt_0 = new TH1D("hPtMatchedJet_hijing_Antikt_0",
					     "For Matching Efficiency, sHIJING, Antikt, R=0.15",
					     NPTBINS, PTBINS);
    hPtMatchedJet_hijing_Antikt_1 = new TH1D("hPtMatchedJet_hijing_Antikt_1",
					     "For Matching Efficiency, sHIJING, Antikt, R=0.2",
					     NPTBINS, PTBINS);
    hPtMatchedJet_hijing_Antikt_2 = new TH1D("hPtMatchedJet_hijing_Antikt_2",
					     "For Matching Efficiency, sHIJING, Antikt, R=0.25",
					     NPTBINS, PTBINS);
    hPtMatchedJet_hijing_Antikt_3 = new TH1D("hPtMatchedJet_hijing_Antikt_3",
					     "For Matching Efficiency, sHIJING, Antikt, R=0.3",
					     NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for (unsigned int i = 0; i < 5; i++)
        {
            hPtTrueVsMatched_hijing_Antikt_0[i] = new TH2D(Form("hPtTrueVsMatched_hijing_Antikt_0_%d", i),
							   "sHIJING, Antikt, R=0.15; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_hijing_Antikt_1[i] = new  TH2D(Form("hPtTrueVsMatched_hijing_Antikt_1_%d", i),
							    "sHIJING, Antikt, R=0.2; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							    NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_hijing_Antikt_2[i] = new  TH2D(Form("hPtTrueVsMatched_hijing_Antikt_2_%d", i),
							    "sHIJING, Antikt, R=0.25; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							    NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_hijing_Antikt_3[i] = new  TH2D(Form("hPtTrueVsMatched_hijing_Antikt_3_%d", i),
							    "sHIJING, Antikt, R=0.3; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							    NPTBINS, PTBINS, NPTBINS, PTBINS);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hForMatchingEfficiency_hijing_Antikt_0[i] = new TH1D(Form("hForMatchingEfficiency_hijing_Antikt_0_%d", i),
								 "For Discriminant Efficiency, sHIJING, Antikt, R=0.15",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_hijing_Antikt_1[i] = new TH1D(Form("hForMatchingEfficiency_hijing_Antikt_1_%d", i),
								 "For Discriminant Efficiency, sHIJING, Antikt, R=0.2",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_hijing_Antikt_2[i] = new TH1D(Form("hForMatchingEfficiency_hijing_Antikt_2_%d", i),
								 "For Discriminant Efficiency, sHIJING, Antikt, R=0.25",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_hijing_Antikt_3[i] = new TH1D(Form("hForMatchingEfficiency_hijing_Antikt_3_%d", i),
								 "For Discriminant Efficiency, sHIJING, Antikt, R=0.3",
								 NPTBINS, PTBINS);

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hPtGoodJet_discriminant_hijing_Antikt_0[i] = new TH1D(Form("hPtGoodJet_discriminant_hijing_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_hijing_Antikt_1[i] = new TH1D(Form("hPtGoodJet_discriminant_hijing_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_hijing_Antikt_2[i] = new TH1D(Form("hPtGoodJet_discriminant_hijing_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_hijing_Antikt_3[i] = new TH1D(Form("hPtGoodJet_discriminant_hijing_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_discriminant_hijing_Antikt_0[i] = new TH1D(Form("hPtMatchedJet_discriminant_hijing_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_hijing_Antikt_1[i] = new TH1D(Form("hPtMatchedJet_discriminant_hijing_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_hijing_Antikt_2[i] = new TH1D(Form("hPtMatchedJet_discriminant_hijing_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_hijing_Antikt_3[i] = new TH1D(Form("hPtMatchedJet_discriminant_hijing_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtGoodJet_centrality_hijing_Antikt_0[i] = new TH1D(Form("hPtGoodJet_centrality_hijing_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_hijing_Antikt_1[i] = new TH1D(Form("hPtGoodJet_centrality_hijing_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_hijing_Antikt_2[i] = new TH1D(Form("hPtGoodJet_centrality_hijing_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_hijing_Antikt_3[i] = new TH1D(Form("hPtGoodJet_centrality_hijing_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_centrality_hijing_Antikt_0[i] = new TH1D(Form("hPtMatchedJet_centrality_hijing_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_hijing_Antikt_1[i] = new TH1D(Form("hPtMatchedJet_centrality_hijing_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_hijing_Antikt_2[i] = new TH1D(Form("hPtMatchedJet_centrality_hijing_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_hijing_Antikt_3[i] = new TH1D(Form("hPtMatchedJet_centrality_hijing_Antikt_3_%d", i), "", NPTBINS, PTBINS);
        }

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Pythia
    //*****************************************************************************************************************************
    //Anti-kt
    //*****************************************************************************************************************************
    hDistance_pythia_Antikt_0 = new TH1D("hDistance_pythia_Antikt_0",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True PYTHIA Jet. Anti-kt, R=0.15",
                                         600, 0.0, 0.6);
    hDistance_pythia_Antikt_1 = new TH1D("hDistance_pythia_Antikt_1",
                                         "#DeltaR between YTHIA+sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True PYTHIA Jet. Anti-kt, R=0.2",
                                         600, 0.0, 0.6);
    hDistance_pythia_Antikt_2 = new TH1D("hDistance_pythia_Antikt_2",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True PYTHIA Jet. Anti-kt, R=0.25",
                                         600, 0.0, 0.6);
    hDistance_pythia_Antikt_3 = new TH1D("hDistance_pythia_Antikt_3",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets (p_{T}>5 && nc>=3) and"
                                         "True PYTHIA Jet. Anti-kt, R=0.3",
                                         600, 0.0, 0.6);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtTrueJet_pythia_Antikt_0 = new TH1D("hPtTrueJet_pythia_Antikt_0",
                                          "p_{T} of PYTHIA True Jet, Antikt, R=0.15; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Antikt_1 = new TH1D("hPtTrueJet_pythia_Antikt_1",
                                          "p_{T} of PYTHIA True Jet, Antikt, R=0.2; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Antikt_2 = new TH1D("hPtTrueJet_pythia_Antikt_2",
                                          "p_{T} of PYTHIA True Jet, Antikt, R=0.25; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Antikt_3 = new TH1D("hPtTrueJet_pythia_Antikt_3",
                                          "p_{T} of PYTHIA True Jet, Antikt, R=0.3; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtMatchedJet_pythia_Antikt_0 = new TH1D("hPtMatchedJet_pythia_Antikt_0",
					     "For Matching Efficiency, PYTHIA+sHIJING, Antikt, R=0.15",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Antikt_1 = new TH1D("hPtMatchedJet_pythia_Antikt_1",
					     "For Matching Efficiency, PYTHIA+sHIJING, Antikt, R=0.2",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Antikt_2 = new TH1D("hPtMatchedJet_pythia_Antikt_2",
					     "For Matching Efficiency, PYTHIA+sHIJING, Antikt, R=0.25",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Antikt_3 = new TH1D("hPtMatchedJet_pythia_Antikt_3",
					     "For Matching Efficiency, PYTHIA+sHIJING, Antikt, R=0.3",
					     NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for (unsigned int i = 0; i < 5; i++)
        {
            hPtTrueVsMatched_pythia_Antikt_0[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Antikt_0_%d", i),
							   "PYTHIA+sHIJING, Antikt, R=0.15; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Antikt_1[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Antikt_1_%d", i),
							   "PYTHIA+sHIJING, Antikt, R=0.2; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Antikt_2[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Antikt_2_%d", i),
							   "PYTHIA+sHIJING, Antikt, R=0.25; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Antikt_3[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Antikt_3_%d", i),
							   "PYTHIA+sHIJING, Antikt, R=0.3; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hForMatchingEfficiency_pythia_Antikt_0[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Antikt_0_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Antikt, R=0.15",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Antikt_1[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Antikt_1_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Antikt, R=0.2",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Antikt_2[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Antikt_2_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Antikt, R=0.25",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Antikt_3[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Antikt_3_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Antikt, R=0.3",
								 NPTBINS, PTBINS);

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hPtGoodJet_discriminant_pythia_Antikt_0[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Antikt_1[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Antikt_2[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Antikt_3[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_discriminant_pythia_Antikt_0[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Antikt_1[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Antikt_2[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Antikt_3[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtGoodJet_centrality_pythia_Antikt_0[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Antikt_1[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Antikt_2[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Antikt_3[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Antikt_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_centrality_pythia_Antikt_0[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Antikt_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Antikt_1[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Antikt_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Antikt_2[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Antikt_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Antikt_3[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Antikt_3_%d", i), "", NPTBINS, PTBINS);
        }

    //*****************************************************************************************************************************
    //Gaussian Filter
    //*****************************************************************************************************************************
    hDistance_pythia_Filter_0 = new TH1D("hDistance_pythia_Filter_0",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets and "
                                         "True PYTHIA Jet. Gaussian Filter, #sigma=0.106",
                                         600, 0.0, 0.6);
    hDistance_pythia_Filter_1 = new TH1D("hDistance_pythia_Filter_1",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets and "
                                         "True PYTHIA Jet. Gaussian Filter, #sigma=0.141",
                                         600, 0.0, 0.6);
    hDistance_pythia_Filter_2 = new TH1D("hDistance_pythia_Filter_2",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets and "
                                         "True PYTHIA Jet. Gausssian Filter, #sigma=0.177",
                                         600, 0.0, 0.6);
    hDistance_pythia_Filter_3 = new TH1D("hDistance_pythia_Filter_3",
                                         "#DeltaR between PYTHIA+sHIJING final particle Jets and "
                                         "True PYTHIA Jet. Gaussian Filter, #sigma=0.212",
                                         600, 0.0, 0.6);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtTrueJet_pythia_Filter_0 = new TH1D("hPtTrueJet_pythia_Filter_0",
                                          "p_{T} of PYTHIA True Jet, Gaussian Filter, #sigma=0.106; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Filter_1 = new TH1D("hPtTrueJet_pythia_Filter_1",
                                          "p_{T} of PYTHIA True Jet, Gaussian Filter, #sigma=0.141; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Filter_2 = new TH1D("hPtTrueJet_pythia_Filter_2",
                                          "p_{T} of PYTHIA True Jet, Gaussian Filter, #sigma=0.177; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    hPtTrueJet_pythia_Filter_3 = new TH1D("hPtTrueJet_pythia_Filter_3",
                                          "p_{T} of PYTHIA True Jet, Gaussian Filter, #sigma=0.212; p_{T, True} (GeV/c)",
                                          NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    hPtMatchedJet_pythia_Filter_0 = new TH1D("hPtMatchedJet_pythia_Filter_0",
					     "For Matching Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.106",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Filter_1 = new TH1D("hPtMatchedJet_pythia_Filter_1",
					     "For Matching Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.141",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Filter_2 = new TH1D("hPtMatchedJet_pythia_Filter_2",
					     "For Matching Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.177",
					     NPTBINS, PTBINS);
    hPtMatchedJet_pythia_Filter_3 = new TH1D("hPtMatchedJet_pythia_Filter_3",
					     "For Matching Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.212",
					     NPTBINS, PTBINS);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for (unsigned int i = 0; i < 5; i++)
        {
            hPtTrueVsMatched_pythia_Filter_0[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Filter_0_%d", i),
							   "PYTHIA+sHIJING, Gaussian Filter, #sigma=0.106; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Filter_1[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Filter_1_%d", i),
							   "PYTHIA+sHIJING, Gaussian Filter, #sigma=0.141; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Filter_2[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Filter_2_%d", i),
							   "PYTHIA+sHIJING, Gaussian Filter, #sigma=0.177; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            hPtTrueVsMatched_pythia_Filter_3[i] = new TH2D(Form("hPtTrueVsMatched_pythia_Filter_3_%d", i),
							   "PYTHIA+sHIJING, Gaussian Filter, #sigma=0.212; p_{T, Matched} (GeV/c);p_{T, True} (GeV/c)",
							   NPTBINS, PTBINS, NPTBINS, PTBINS);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hForMatchingEfficiency_pythia_Filter_0[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Filter_0_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.106",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Filter_1[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Filter_1_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.141",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Filter_2[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Filter_2_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.177",
								 NPTBINS, PTBINS);
            hForMatchingEfficiency_pythia_Filter_3[i] = new TH1D(Form("hForMatchingEfficiency_pythia_Filter_3_%d", i),
								 "For Discriminant Efficiency, PYTHIA+sHIJING, Gaussian Filter, #sigma=0.212",
								 NPTBINS, PTBINS);

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hPtGoodJet_discriminant_pythia_Filter_0[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Filter_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Filter_1[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Filter_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Filter_2[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Filter_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_discriminant_pythia_Filter_3[i] = new TH1D(Form("hPtGoodJet_discriminant_pythia_Filter_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_discriminant_pythia_Filter_0[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Filter_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Filter_1[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Filter_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Filter_2[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Filter_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_discriminant_pythia_Filter_3[i] = new TH1D(Form("hPtMatchedJet_discriminant_pythia_Filter_3_%d", i), "", NPTBINS, PTBINS);

            hPtGoodJet_centrality_pythia_Filter_0[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Filter_0_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Filter_1[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Filter_1_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Filter_2[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Filter_2_%d", i), "", NPTBINS, PTBINS);
            hPtGoodJet_centrality_pythia_Filter_3[i] = new TH1D(Form("hPtGoodJet_centrality_pythia_Filter_3_%d", i), "", NPTBINS, PTBINS);

            hPtMatchedJet_centrality_pythia_Filter_0[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Filter_0_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Filter_1[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Filter_1_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Filter_2[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Filter_2_%d", i), "", NPTBINS, PTBINS);
            hPtMatchedJet_centrality_pythia_Filter_3[i] = new TH1D(Form("hPtMatchedJet_centrality_pythia_Filter_3_%d", i), "", NPTBINS, PTBINS);

        }

    if(fillTrees)
        {
            //*****************************************************************************************************************************
            //*****************************************************************************************************************************
            //nTuples
            //*****************************************************************************************************************************
            //sHijing
            //*****************************************************************************************************************************
            finalHijing = new TTree("finalHijing", "Final particles of sHIJING");
            finalHijing->Branch("finalHijing", &hijing_particle_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_true_antikt_jet_0 = new TTree("hijing_true_antikt_jet_0", "True Jets of sHIJING, Anti-kt, R=0.15");
            hijing_true_antikt_jet_0->Branch("hijing_true_antikt_jet_0", &hijing_true_antikt_jet_0_list);
            hijing_true_antikt_jet_1 = new TTree("hijing_true_antikt_jet_1", "True Jets of sHIJING, Anti-kt, R=0.2");
            hijing_true_antikt_jet_1->Branch("hijing_true_antikt_jet_1", &hijing_true_antikt_jet_1_list);
            hijing_true_antikt_jet_2 = new TTree("hijing_true_antikt_jet_2", "True Jets of sHIJING, Anti-kt, R=0.25");
            hijing_true_antikt_jet_2->Branch("hijing_true_antikt_jet_2", &hijing_true_antikt_jet_2_list);
            hijing_true_antikt_jet_3 = new TTree("hijing_true_antikt_jet_3", "True Jets of sHIJING, Anti-kt, R=0.3");
            hijing_true_antikt_jet_3->Branch("hijing_true_antikt_jet_3", &hijing_true_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_good_antikt_jet_0 = new TTree("hijing_good_antikt_jet_0",
                                                 "Reconstructed Jets (p_{T}>5 && nc>=3) from all sHIJING particles, anti-kt, R = 0.15");
            hijing_good_antikt_jet_0->Branch("hijing_good_antikt_jet_0", &hijing_good_antikt_jet_0_list);
            hijing_good_antikt_jet_1 = new TTree("hijing_good_antikt_jet_1",
                                                 "Reconstructed Jets (p_{T}>5 && nc>=3) from all sHIJING particles, anti-kt, R = 0.2");
            hijing_good_antikt_jet_1->Branch("hijing_good_antikt_jet_1", &hijing_good_antikt_jet_1_list);
            hijing_good_antikt_jet_2 = new TTree("hijing_good_antikt_jet_2",
                                                 "Reconstructed Jets (p_{T}>5 && nc>=3) from all sHIJING particles, anti-kt, R = 0.25");
            hijing_good_antikt_jet_2->Branch("hijing_good_antikt_jet_2", &hijing_good_antikt_jet_2_list);
            hijing_good_antikt_jet_3 = new TTree("hijing_good_antikt_jet_3",
                                                 "Reconstructed Jets (p_{T}>5 && nc>=3) from all sHIJING particles, anti-kt, R = 0.3");
            hijing_good_antikt_jet_3->Branch("hijing_good_antikt_jet_3", &hijing_good_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_good_antikt_jet_constituents_0 = new TTree("hijing_good_antikt_jet_constituents_0",
							      "Constituents of sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.15");
            hijing_good_antikt_jet_constituents_0->Branch("hijing_good_antikt_jet_constituents_0", &hijing_good_antikt_jet_constituents_0_list);
            hijing_good_antikt_jet_constituents_1 = new TTree("hijing_good_antikt_jet_constituents_1",
							      "Constituents of sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.2");
            hijing_good_antikt_jet_constituents_1->Branch("hijing_good_antikt_jet_constituents_1", &hijing_good_antikt_jet_constituents_1_list);
            hijing_good_antikt_jet_constituents_2 = new TTree("hijing_good_antikt_jet_constituents_2",
							      "Constituents of sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.25");
            hijing_good_antikt_jet_constituents_2->Branch("hijing_good_antikt_jet_constituents_2", &hijing_good_antikt_jet_constituents_2_list);
            hijing_good_antikt_jet_constituents_3 = new TTree("hijing_good_antikt_jet_constituents_3",
							      "Constituents of sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.3");
            hijing_good_antikt_jet_constituents_3->Branch("hijing_good_antikt_jet_constituents_3", &hijing_good_antikt_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_matched_antikt_jet_0 = new TTree("hijing_matched_antikt_jet_0",
                                                    "Matched Jets of sHIJING, anti-kt, R = 0.15");
            hijing_matched_antikt_jet_0->Branch("hijing_matched_antikt_jet_0", &hijing_matched_antikt_jet_0_list);
            hijing_matched_antikt_jet_1 = new TTree("hijing_matched_antikt_jet_1",
                                                    "Matched Jets of sHIJING, anti-kt, R = 0.2");
            hijing_matched_antikt_jet_1->Branch("hijing_matched_antikt_jet_1", &hijing_matched_antikt_jet_1_list);
            hijing_matched_antikt_jet_2 = new TTree("hijing_matched_antikt_jet_2",
                                                    "Matched Jets of sHIJING, anti-kt, R = 0.25");
            hijing_matched_antikt_jet_2->Branch("hijing_matched_antikt_jet_2", &hijing_matched_antikt_jet_2_list);
            hijing_matched_antikt_jet_3 = new TTree("hijing_matched_antikt_jet_3",
                                                    "Matched Jets of sHIJING, anti-kt, R = 0.3");
            hijing_matched_antikt_jet_3->Branch("hijing_matched_antikt_jet_3", &hijing_matched_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_matched_antikt_jet_constituents_0 = new TTree("hijing_matched_jet_constituents_0",
								 "Constituents of Matched sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.15");
            hijing_matched_antikt_jet_constituents_0->Branch("hijing_matched_jet_constituents_0", &hijing_matched_antikt_jet_constituents_0_list);
            hijing_matched_antikt_jet_constituents_1 = new TTree("hijing_matched_jet_constituents_1",
								 "Constituents of Matched sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.2");
            hijing_matched_antikt_jet_constituents_1->Branch("hijing_matched_jet_constituents_1", &hijing_matched_antikt_jet_constituents_1_list);
            hijing_matched_antikt_jet_constituents_2 = new TTree("hijing_matched_jet_constituents_2",
								 "Constituents of Matched sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.25");
            hijing_matched_antikt_jet_constituents_2->Branch("hijing_matched_jet_constituents_2", &hijing_matched_antikt_jet_constituents_2_list);
            hijing_matched_antikt_jet_constituents_3 = new TTree("hijing_matched_jet_constituents_3",
								 "Constituents of Matched sHIJING Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.3");
            hijing_matched_antikt_jet_constituents_3->Branch("hijing_matched_jet_constituents_3", &hijing_matched_antikt_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            hijing_matched_jet_subtracted_0 = new TTree("hijing_matched_jet_subtracted_0",
							"Constituents of Matched sHIJING Jet constituents subtraced list, Anti-kt, R=0.15");
            hijing_matched_jet_subtracted_0->Branch("hijing_matched_jet_subtracted_0", &hijing_matched_jet_subtracted_0_list);
            hijing_matched_jet_subtracted_1 = new TTree("hijing_matched_jet_subtracted_1",
							"Constituents of Matched sHIJING Jet constituents subtraced list, Anti-kt, R=0.2");
            hijing_matched_jet_subtracted_1->Branch("hijing_matched_jet_subtracted_1", &hijing_matched_jet_subtracted_1_list);
            hijing_matched_jet_subtracted_2 = new TTree("hijing_matched_jet_subtracted_2",
							"Constituents of Matched sHIJING Jet constituents subtraced list, Anti-kt, R=0.25");
            hijing_matched_jet_subtracted_2->Branch("hijing_matched_jet_subtracted_2", &hijing_matched_jet_subtracted_2_list);
            hijing_matched_jet_subtracted_3 = new TTree("hijing_matched_jet_subtracted_3",
							"Constituents of Matched sHIJING Jets constituents subtraced list, Anti-kt, R=0.3");
            hijing_matched_jet_subtracted_3->Branch("hijing_matched_jet_subtracted_3", &hijing_matched_jet_subtracted_3_list);

            //*****************************************************************************************************************************
            //Pythia
            finalPythia = new TTree("finalPythia", "Final particles of Pythia");
            finalPythia->Branch("finalPythia", &pythia_particle_list);
            //*****************************************************************************************************************************

            //*****************************************************************************************************************************
            //Anti-kt
            //*****************************************************************************************************************************
            pythia_true_antikt_jet_0 = new TTree("pythia_true_antikt_jet_0", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Anti-kt, R=0.15");
            pythia_true_antikt_jet_0->Branch("pythia_true_antikt_jet_0", &pythia_true_antikt_jet_0_list);
            pythia_true_antikt_jet_1 = new TTree("pythia_true_antikt_jet_1", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Anti-kt, R=0.2");
            pythia_true_antikt_jet_1->Branch("pythia_true_antikt_jet_1", &pythia_true_antikt_jet_1_list);
            pythia_true_antikt_jet_2 = new TTree("pythia_true_antikt_jet_2", "True Jets of PYTHIA, Anti-kt, R=0.25");
            pythia_true_antikt_jet_2->Branch("pythia_true_antikt_jet_2", &pythia_true_antikt_jet_2_list);
            pythia_true_antikt_jet_3 = new TTree("pythia_true_antikt_jet_3", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Anti-kt, R=0.3");
            pythia_true_antikt_jet_3->Branch("pythia_true_antikt_jet_3", &pythia_true_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            pythia_true_antikt_jet_constituents_0 = new TTree("pythia_true_antikt_jet_constituents_0",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.15");
            pythia_true_antikt_jet_constituents_0->Branch("pythia_true_antikt_jet_constituents_0", &pythia_true_antikt_jet_constituents_0_list);
            pythia_true_antikt_jet_constituents_1 = new TTree("pythia_true_antikt_jet_constituents_1",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.2");
            pythia_true_antikt_jet_constituents_1->Branch("pythia_true_antikt_jet_constituents_1", &pythia_true_antikt_jet_constituents_1_list);
            pythia_true_antikt_jet_constituents_2 = new TTree("pythia_true_antikt_jet_constituents_2",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.25");
            pythia_true_antikt_jet_constituents_2->Branch("pythia_true_antikt_jet_constituents_2", &pythia_true_antikt_jet_constituents_2_list);
            pythia_true_antikt_jet_constituents_3 = new TTree("pythia_true_antikt_jet_constituents_3",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Anti-kt, R=0.3");
            pythia_true_antikt_jet_constituents_3->Branch("pythia_true_antikt_jet_constituents_3", &pythia_true_antikt_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_particle_antikt_0 = new TTree("embedded_particle_antikt_0", "PYTHIA+sHIJING, Anti-kt, R = 0.15");
            embedded_particle_antikt_0->Branch("embedded_particle_antikt_0", &embedded_particle_antikt_0_list);
            embedded_particle_antikt_1 = new TTree("embedded_particle_antikt_1", "PYTHIA+sHIJING, Anti-kt, R = 0.2");
            embedded_particle_antikt_1->Branch("embedded_particle_antikt_1", &embedded_particle_antikt_1_list);
            embedded_particle_antikt_2 = new TTree("embedded_particle_antikt_2", "PYTHIA+sHIJING, Anti-kt, R = 0.25");
            embedded_particle_antikt_2->Branch("embedded_particle_antikt_2", &embedded_particle_antikt_2_list);
            embedded_particle_antikt_3 = new TTree("embedded_particle_antikt_3", "PYTHIA+sHIJING, Anti-kt, R = 0.3");
            embedded_particle_antikt_3->Branch("embedded_particle_antikt_3", &embedded_particle_antikt_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_good_antikt_jet_0 = new TTree("embedded_good_antikt_jet_0", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.15");
            embedded_good_antikt_jet_0->Branch("embedded_good_antikt_jet_0", &embedded_good_antikt_jet_0_list);
            embedded_good_antikt_jet_1 = new TTree("embedded_good_antikt_jet_1", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.2");
            embedded_good_antikt_jet_1->Branch("embedded_good_antikt_jet_1", &embedded_good_antikt_jet_1_list);
            embedded_good_antikt_jet_2 = new TTree("embedded_good_antikt_jet_2", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.25");
            embedded_good_antikt_jet_2->Branch("embedded_good_antikt_jet_2", &embedded_good_antikt_jet_2_list);
            embedded_good_antikt_jet_3 = new TTree("embedded_good_antikt_jet_3", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.3");
            embedded_good_antikt_jet_3->Branch("embedded_good_antikt_jet_3", &embedded_good_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_good_antikt_jet_constituents_0 = new TTree("embedded_good_antikt_jet_constituents_0",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.15");
            embedded_good_antikt_jet_constituents_0->Branch("embedded_good_antikt_jet_constituents_0", &embedded_good_antikt_jet_constituents_0_list);
            embedded_good_antikt_jet_constituents_1 = new TTree("embedded_good_antikt_jet_constituents_1",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.2");
            embedded_good_antikt_jet_constituents_1->Branch("embedded_good_antikt_jet_constituents_1", &embedded_good_antikt_jet_constituents_1_list);
            embedded_good_antikt_jet_constituents_2 = new TTree("embedded_good_antikt_jet_constituents_2",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.25");
            embedded_good_antikt_jet_constituents_2->Branch("embedded_good_antikt_jet_constituents_2", &embedded_good_antikt_jet_constituents_2_list);
            embedded_good_antikt_jet_constituents_3 = new TTree("embedded_good_antikt_jet_constituents_3",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Anti-kt, R = 0.3");
            embedded_good_antikt_jet_constituents_3->Branch("embedded_good_antikt_jet_constituents_3", &embedded_good_antikt_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_matched_antikt_jet_0 = new TTree("embedded_matched_antikt_jet_0",
						      "Matched Jets of PYTHIA+sHIJING, anti-kt, R = 0.15");
            embedded_matched_antikt_jet_0->Branch("embedded_matched_antikt_jet_0", &embedded_matched_antikt_jet_0_list);
            embedded_matched_antikt_jet_1 = new TTree("embedded_matched_antikt_jet_1",
						      "Matched Jets of PYTHIA+sHIJING, anti-kt, R = 0.2");
            embedded_matched_antikt_jet_1->Branch("embedded_matched_antikt_jet_1", &embedded_matched_antikt_jet_1_list);
            embedded_matched_antikt_jet_2 = new TTree("embedded_matched_antikt_jet_2",
						      "Matched Jets of PYTHIA+sHIJING, anti-kt, R = 0.25");
            embedded_matched_antikt_jet_2->Branch("embedded_matched_antikt_jet_2", &embedded_matched_antikt_jet_2_list);
            embedded_matched_antikt_jet_3 = new TTree("embedded_matched_antikt_jet_3",
						      "Matched Jets of PYTHIA+sHIJING, anti-kt, R = 0.3");
            embedded_matched_antikt_jet_3->Branch("embedded_matched_antikt_jet_3", &embedded_matched_antikt_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            embedded_matched_antikt_jet_constituents_0 = new TTree("embedded_matched_antikt_jet_constituents_0",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Anti-kt, R=0.15");
            embedded_matched_antikt_jet_constituents_0->Branch("embedded_matched_antikt_jet_constituents_0", &embedded_matched_antikt_jet_constituents_0_list);
            embedded_matched_antikt_jet_constituents_1 = new TTree("embedded_matched_antikt_jet_constituents_1",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Anti-kt, R=0.2");
            embedded_matched_antikt_jet_constituents_1->Branch("embedded_matched_antikt_jet_constituents_1", &embedded_matched_antikt_jet_constituents_1_list);
            embedded_matched_antikt_jet_constituents_2 = new TTree("embedded_matched_antikt_jet_constituents_2",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Anti-kt, R=0.25");
            embedded_matched_antikt_jet_constituents_2->Branch("embedded_matched_antikt_jet_constituents_2", &embedded_matched_antikt_jet_constituents_2_list);
            embedded_matched_antikt_jet_constituents_3 = new TTree("embedded_matched_antikt_jet_constituents_3",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Anti-kt, R=0.3");
            embedded_matched_antikt_jet_constituents_3->Branch("embedded_matched_antikt_jet_constituents_3", &embedded_matched_antikt_jet_constituents_3_list);
            //*****************************************************************************************************************************
            //Gaussian Filter
            //*****************************************************************************************************************************
            pythia_true_filter_jet_0 = new TTree("pythia_true_filter_jet_0", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Gaussian Filter, #sigma=0.106");
            pythia_true_filter_jet_0->Branch("pythia_true_filter_jet_0", &pythia_true_filter_jet_0_list);
            pythia_true_filter_jet_1 = new TTree("pythia_true_filter_jet_1", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Gaussian Filter, #sigma=0.141");
            pythia_true_filter_jet_1->Branch("pythia_true_filter_jet_1", &pythia_true_filter_jet_1_list);
            pythia_true_filter_jet_2 = new TTree("pythia_true_filter_jet_2", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Gaussian Filter, #sigma=0.177");
            pythia_true_filter_jet_2->Branch("pythia_true_filter_jet_2", &pythia_true_filter_jet_2_list);
            pythia_true_filter_jet_3 = new TTree("pythia_true_filter_jet_3", "True Jets (p_{T}>5 && nc>=3) of PYTHIA, Gaussian Filter, #sigma=0.212");
            pythia_true_filter_jet_3->Branch("pythia_true_filter_jet_3", &pythia_true_filter_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            pythia_true_filter_jet_constituents_0 = new TTree("pythia_true_filter_jet_constituents_0",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Gaussian Filter, #sigma=0.106");
            pythia_true_filter_jet_constituents_0->Branch("pythia_true_filter_jet_constituents_0", &pythia_true_filter_jet_constituents_0_list);
            pythia_true_filter_jet_constituents_1 = new TTree("pythia_true_filter_jet_constituents_1",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Gaussian Filter, #sigma=0.141");
            pythia_true_filter_jet_constituents_1->Branch("pythia_true_filter_jet_constituents_1", &pythia_true_filter_jet_constituents_1_list);
            pythia_true_filter_jet_constituents_2 = new TTree("pythia_true_filter_jet_constituents_2",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Gaussian Filter, #sigma=0.177");
            pythia_true_filter_jet_constituents_2->Branch("pythia_true_filter_jet_constituents_2", &pythia_true_filter_jet_constituents_2_list);
            pythia_true_filter_jet_constituents_3 = new TTree("pythia_true_filter_jet_constituents_3",
							      "Constituents of True PYTHIA Jets (p_{T}>5 && nc>=3), Gaussian Filter, #sigma=0.212");
            pythia_true_filter_jet_constituents_3->Branch("pythia_true_filter_jet_constituents_3", &pythia_true_filter_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_particle_filter_0 = new TTree("embedded_particle_filter_0", "PYTHIA+sHIJING, Gaussian Filter, #sigma = 0.106");
            embedded_particle_filter_0->Branch("embedded_particle_filter_0", &embedded_particle_filter_0_list);
            embedded_particle_filter_1 = new TTree("embedded_particle_filter_1", "PYTHIA+sHIJING, Gaussian Filter, #sigma = 0.141");
            embedded_particle_filter_1->Branch("embedded_particle_filter_1", &embedded_particle_filter_1_list);
            embedded_particle_filter_2 = new TTree("embedded_particle_filter_2", "PYTHIA+sHIJING, Gaussian Filter, #sigma = 0.177");
            embedded_particle_filter_2->Branch("embedded_particle_filter_2", &embedded_particle_filter_2_list);
            embedded_particle_filter_3 = new TTree("embedded_particle_filter_3", "PYTHIA+sHIJING, Gaussian Filter, #sigma = 0.212");
            embedded_particle_filter_3->Branch("embedded_particle_filter_3", &embedded_particle_filter_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_good_filter_jet_0 = new TTree("embedded_good_filter_jet_0", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.106");
            embedded_good_filter_jet_0->Branch("embedded_good_filter_jet_0", &embedded_good_filter_jet_0_list);
            embedded_good_filter_jet_1 = new TTree("embedded_good_filter_jet_1", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.141");
            embedded_good_filter_jet_1->Branch("embedded_good_filter_jet_1", &embedded_good_filter_jet_1_list);
            embedded_good_filter_jet_2 = new TTree("embedded_good_filter_jet_2", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.177");
            embedded_good_filter_jet_2->Branch("embedded_good_filter_jet_2", &embedded_good_filter_jet_2_list);
            embedded_good_filter_jet_3 = new TTree("embedded_good_filter_jet_3", "PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.212");
            embedded_good_filter_jet_3->Branch("embedded_good_filter_jet_3", &embedded_good_filter_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_good_filter_jet_constituents_0 = new TTree("embedded_good_filter_jet_constituents_0",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.106");
            embedded_good_filter_jet_constituents_0->Branch("embedded_good_filter_jet_constituents_0", &embedded_good_filter_jet_constituents_0_list);
            embedded_good_filter_jet_constituents_1 = new TTree("embedded_good_filter_jet_constituents_1",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.141");
            embedded_good_filter_jet_constituents_1->Branch("embedded_good_filter_jet_constituents_1", &embedded_good_filter_jet_constituents_1_list);
            embedded_good_filter_jet_constituents_2 = new TTree("embedded_good_filter_jet_constituents_2",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.177");
            embedded_good_filter_jet_constituents_2->Branch("embedded_good_filter_jet_constituents_2", &embedded_good_filter_jet_constituents_2_list);
            embedded_good_filter_jet_constituents_3 = new TTree("embedded_good_filter_jet_constituents_3",
								"Constituents of PYTHIA+sHIJING Jet (p_{T}>5 && nc>=3), Gaussian Filter, #sigma = 0.212");
            embedded_good_filter_jet_constituents_3->Branch("embedded_good_filter_jet_constituents_3", &embedded_good_filter_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_matched_filter_jet_0 = new TTree("embedded_matched_filter_jet_0",
						      "Matched Jets of PYTHIA+sHIJING, Gaussian Filter, #sigma=0.106");
            embedded_matched_filter_jet_0->Branch("embedded_matched_filter_jet_0", &embedded_matched_filter_jet_0_list);
            embedded_matched_filter_jet_1 = new TTree("embedded_matched_filter_jet_1",
						      "Matched Jets of PYTHIA+sHIJING, Gaussian Filter, #sigma=0.141");
            embedded_matched_filter_jet_1->Branch("embedded_matched_filter_jet_1", &embedded_matched_filter_jet_1_list);
            embedded_matched_filter_jet_2 = new TTree("embedded_matched_filter_jet_2",
						      "Matched Jets of PYTHIA+sHIJING, Gaussian Filter, #sigma=0.177");
            embedded_matched_filter_jet_2->Branch("embedded_matched_filter_jet_2", &embedded_matched_filter_jet_2_list);
            embedded_matched_filter_jet_3 = new TTree("embedded_matched_filter_jet_3",
						      "Matched Jets of PYTHIA+sHIJING, Gaussian Filter, #sigma=0.212");
            embedded_matched_filter_jet_3->Branch("embedded_matched_filter_jet_3", &embedded_matched_filter_jet_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            embedded_matched_filter_jet_constituents_0 = new TTree("embedded_matched_filter_jet_constituents_0",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Gaussian Filter, #sigma=0.106");
            embedded_matched_filter_jet_constituents_0->Branch("embedded_matched_filter_jet_constituents_0", &embedded_matched_filter_jet_constituents_0_list);
            embedded_matched_filter_jet_constituents_1 = new TTree("embedded_matched_filter_jet_constituents_1",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Gaussian Filter, #sigma=0.141");
            embedded_matched_filter_jet_constituents_1->Branch("embedded_matched_filter_jet_constituents_1", &embedded_matched_filter_jet_constituents_1_list);
            embedded_matched_filter_jet_constituents_2 = new TTree("embedded_matched_filter_jet_constituents_2",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Gaussian Filter, #sigma=0.177");
            embedded_matched_filter_jet_constituents_2->Branch("embedded_matched_filter_jet_constituents_2", &embedded_matched_filter_jet_constituents_2_list);
            embedded_matched_filter_jet_constituents_3 = new TTree("embedded_matched_filter_jet_constituents_3",
								   "Constituents of Matched PYTHIA+sHIJING Jets, Gaussian Filter, #sigma=0.212");
            embedded_matched_filter_jet_constituents_3->Branch("embedded_matched_filter_jet_constituents_3", &embedded_matched_filter_jet_constituents_3_list);
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }

    return 0;
}

int JetSimWithoutDetector::InitRun(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }
    nRunEvents = 0;

    return EVENT_OK;
}

int JetSimWithoutDetector::ResetEvent(PHCompositeNode *topNode)
{
    //*****************************************************************************************************************************
    //sHijing
    //*****************************************************************************************************************************
    hijing_particle_list.clear();

    hijing_true_antikt_jet_0_list.clear();
    hijing_true_antikt_jet_1_list.clear();
    hijing_true_antikt_jet_2_list.clear();
    hijing_true_antikt_jet_3_list.clear();

    hijing_good_antikt_jet_0_list.clear();
    hijing_good_antikt_jet_1_list.clear();
    hijing_good_antikt_jet_2_list.clear();
    hijing_good_antikt_jet_3_list.clear();

    hijing_good_antikt_jet_constituents_0_list.clear();
    hijing_good_antikt_jet_constituents_1_list.clear();
    hijing_good_antikt_jet_constituents_2_list.clear();
    hijing_good_antikt_jet_constituents_3_list.clear();

    hijing_matched_antikt_jet_0_list.clear();
    hijing_matched_antikt_jet_1_list.clear();
    hijing_matched_antikt_jet_2_list.clear();
    hijing_matched_antikt_jet_3_list.clear();

    hijing_matched_antikt_jet_constituents_0_list.clear();
    hijing_matched_antikt_jet_constituents_1_list.clear();
    hijing_matched_antikt_jet_constituents_2_list.clear();
    hijing_matched_antikt_jet_constituents_3_list.clear();

    hijing_matched_jet_subtracted_0_list.clear();
    hijing_matched_jet_subtracted_1_list.clear();
    hijing_matched_jet_subtracted_2_list.clear();
    hijing_matched_jet_subtracted_3_list.clear();

    //*****************************************************************************************************************************
    //Pythia
    pythia_particle_list.clear();
    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Anti-kt
    //*****************************************************************************************************************************
    pythia_true_antikt_jet_0_list.clear();
    pythia_true_antikt_jet_1_list.clear();
    pythia_true_antikt_jet_2_list.clear();
    pythia_true_antikt_jet_3_list.clear();

    pythia_true_antikt_jet_constituents_0_list.clear();
    pythia_true_antikt_jet_constituents_1_list.clear();
    pythia_true_antikt_jet_constituents_2_list.clear();
    pythia_true_antikt_jet_constituents_3_list.clear();

    embedded_particle_antikt_0_list.clear();
    embedded_particle_antikt_1_list.clear();
    embedded_particle_antikt_2_list.clear();
    embedded_particle_antikt_3_list.clear();

    embedded_good_antikt_jet_0_list.clear();
    embedded_good_antikt_jet_1_list.clear();
    embedded_good_antikt_jet_2_list.clear();
    embedded_good_antikt_jet_3_list.clear();

    embedded_good_antikt_jet_constituents_0_list.clear();
    embedded_good_antikt_jet_constituents_1_list.clear();
    embedded_good_antikt_jet_constituents_2_list.clear();
    embedded_good_antikt_jet_constituents_3_list.clear();

    embedded_matched_antikt_jet_0_list.clear();
    embedded_matched_antikt_jet_1_list.clear();
    embedded_matched_antikt_jet_2_list.clear();
    embedded_matched_antikt_jet_3_list.clear();

    embedded_matched_antikt_jet_constituents_0_list.clear();
    embedded_matched_antikt_jet_constituents_1_list.clear();
    embedded_matched_antikt_jet_constituents_2_list.clear();
    embedded_matched_antikt_jet_constituents_3_list.clear();

    //*****************************************************************************************************************************
    //Gaussian Filter
    //*****************************************************************************************************************************
    pythia_true_filter_jet_0_list.clear();
    pythia_true_filter_jet_1_list.clear();
    pythia_true_filter_jet_2_list.clear();
    pythia_true_filter_jet_3_list.clear();

    pythia_true_filter_jet_constituents_0_list.clear();
    pythia_true_filter_jet_constituents_1_list.clear();
    pythia_true_filter_jet_constituents_2_list.clear();
    pythia_true_filter_jet_constituents_3_list.clear();

    embedded_particle_filter_0_list.clear();
    embedded_particle_filter_1_list.clear();
    embedded_particle_filter_2_list.clear();
    embedded_particle_filter_3_list.clear();

    embedded_good_filter_jet_0_list.clear();
    embedded_good_filter_jet_1_list.clear();
    embedded_good_filter_jet_2_list.clear();
    embedded_good_filter_jet_3_list.clear();

    embedded_good_filter_jet_constituents_0_list.clear();
    embedded_good_filter_jet_constituents_1_list.clear();
    embedded_good_filter_jet_constituents_2_list.clear();
    embedded_good_filter_jet_constituents_3_list.clear();

    embedded_matched_filter_jet_0_list.clear();
    embedded_matched_filter_jet_1_list.clear();
    embedded_matched_filter_jet_2_list.clear();
    embedded_matched_filter_jet_3_list.clear();

    embedded_matched_filter_jet_constituents_0_list.clear();
    embedded_matched_filter_jet_constituents_1_list.clear();
    embedded_matched_filter_jet_constituents_2_list.clear();
    embedded_matched_filter_jet_constituents_3_list.clear();


    return EVENT_OK;
}


int JetSimWithoutDetector::process_event(PHCompositeNode *topNode)
{
    nHijingMultiplicity = 0;
    // Get the data I need...
    // Get PYTHIA Particles
    phpythia                              = getClass<PHPythiaContainer>          (topNode, "PHPythia");
    if (!phpythia)
        {
            cout << "No PHPythia! No sense continuing" << endl;
            exit(1);
        }

    //Find the HepMC data
    PHHepMCGenEvent *genEvt = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
    shijing =  genEvt->getEvent();
    if (!shijing)
        {
            cout << "No sHijing! No sense continuing" << endl;
            exit(1);
        }

    // Informational message...
    nRunEvents++;
    if (nRunEvents % 1000 == 0 && verbosity)
        {
            if (verbo)
                {
                    cout << "Events for run " << runNumber << " = " << nRunEvents << endl;
                }
        }

    //Global things
    HepMC::HeavyIon* heavyIon = shijing->heavy_ion();
    impactParameter = (float)heavyIon->impact_parameter();
    eventPlane = (float)heavyIon->event_plane_angle();

    hImpactParameter->Fill(impactParameter);
    hEventPlane->Fill(eventPlane);

    //*****************************************************************************************************************************
    //sHijing
    //*****************************************************************************************************************************
    GetsHijing();

    //Determine centrality by Total Final State Particles
    hForCentrality->Fill(nHijingMultiplicity);
    //centralityBin = getCentralityBin(nHijingMultiplicity);

    //Get centrality from Impact Parameter
    centralityBin = getCentralityBin(impactParameter);

    GetAntiKt(hijing_particle_list, hijing_good_antikt_jet_0_list, hijing_good_antikt_jet_constituents_0_list, antikt_00);
    GetAntiKt(hijing_particle_list, hijing_good_antikt_jet_1_list, hijing_good_antikt_jet_constituents_1_list, antikt_01);
    GetAntiKt(hijing_particle_list, hijing_good_antikt_jet_2_list, hijing_good_antikt_jet_constituents_2_list, antikt_02);
    GetAntiKt(hijing_particle_list, hijing_good_antikt_jet_3_list, hijing_good_antikt_jet_constituents_3_list, antikt_03);

    GetMatchedJets(hijing_true_antikt_jet_0_list, hijing_good_antikt_jet_0_list, hijing_good_antikt_jet_constituents_0_list,
                   hijing_matched_antikt_jet_0_list, hijing_matched_antikt_jet_constituents_0_list,
                   hDistance_hijing_Antikt_0, hPtTrueJet_hijing_Antikt_0, hPtTrueVsMatched_hijing_Antikt_0,
                   hPtMatchedJet_hijing_Antikt_0, hForMatchingEfficiency_hijing_Antikt_0,
                   minDeltaRHijing0);
    GetMatchedJets(hijing_true_antikt_jet_1_list, hijing_good_antikt_jet_1_list, hijing_good_antikt_jet_constituents_1_list,
                   hijing_matched_antikt_jet_1_list, hijing_matched_antikt_jet_constituents_1_list,
                   hDistance_hijing_Antikt_1, hPtTrueJet_hijing_Antikt_1, hPtTrueVsMatched_hijing_Antikt_1,
                   hPtMatchedJet_hijing_Antikt_1, hForMatchingEfficiency_hijing_Antikt_1,
                   minDeltaRHijing1);
    GetMatchedJets(hijing_true_antikt_jet_2_list, hijing_good_antikt_jet_2_list, hijing_good_antikt_jet_constituents_2_list,
                   hijing_matched_antikt_jet_2_list, hijing_matched_antikt_jet_constituents_2_list,
                   hDistance_hijing_Antikt_2, hPtTrueJet_hijing_Antikt_2, hPtTrueVsMatched_hijing_Antikt_2,
                   hPtMatchedJet_hijing_Antikt_2, hForMatchingEfficiency_hijing_Antikt_2,
                   minDeltaRHijing2);
    GetMatchedJets(hijing_true_antikt_jet_3_list, hijing_good_antikt_jet_3_list, hijing_good_antikt_jet_constituents_3_list,
                   hijing_matched_antikt_jet_3_list, hijing_matched_antikt_jet_constituents_3_list,
                   hDistance_hijing_Antikt_3, hPtTrueJet_hijing_Antikt_3, hPtTrueVsMatched_hijing_Antikt_3,
                   hPtMatchedJet_hijing_Antikt_3, hForMatchingEfficiency_hijing_Antikt_3,
                   minDeltaRHijing3);

    SubtractMatchedJetConstituents(hijing_particle_list, hijing_matched_antikt_jet_constituents_0_list, hijing_matched_jet_subtracted_0_list);
    SubtractMatchedJetConstituents(hijing_particle_list, hijing_matched_antikt_jet_constituents_1_list, hijing_matched_jet_subtracted_1_list);
    SubtractMatchedJetConstituents(hijing_particle_list, hijing_matched_antikt_jet_constituents_2_list, hijing_matched_jet_subtracted_2_list);
    SubtractMatchedJetConstituents(hijing_particle_list, hijing_matched_antikt_jet_constituents_3_list, hijing_matched_jet_subtracted_3_list);

    FillHistograms(hijing_good_antikt_jet_0_list, hijing_matched_antikt_jet_0_list,
                   hPtGoodJet_discriminant_hijing_Antikt_0, hPtMatchedJet_discriminant_hijing_Antikt_0,
                   hPtGoodJet_centrality_hijing_Antikt_0, hPtMatchedJet_centrality_hijing_Antikt_0);
    FillHistograms(hijing_good_antikt_jet_1_list, hijing_matched_antikt_jet_1_list,
                   hPtGoodJet_discriminant_hijing_Antikt_1, hPtMatchedJet_discriminant_hijing_Antikt_1,
                   hPtGoodJet_centrality_hijing_Antikt_1, hPtMatchedJet_centrality_hijing_Antikt_1);
    FillHistograms(hijing_good_antikt_jet_2_list, hijing_matched_antikt_jet_2_list,
                   hPtGoodJet_discriminant_hijing_Antikt_2, hPtMatchedJet_discriminant_hijing_Antikt_2,
                   hPtGoodJet_centrality_hijing_Antikt_2, hPtMatchedJet_centrality_hijing_Antikt_2);
    FillHistograms(hijing_good_antikt_jet_3_list, hijing_matched_antikt_jet_3_list,
                   hPtGoodJet_discriminant_hijing_Antikt_3, hPtMatchedJet_discriminant_hijing_Antikt_3,
                   hPtGoodJet_centrality_hijing_Antikt_3, hPtMatchedJet_centrality_hijing_Antikt_3);

    //*****************************************************************************************************************************
    //Pythia
    //*****************************************************************************************************************************
    GetPythia();
    //*****************************************************************************************************************************
    //Anti-kt
    //*****************************************************************************************************************************
    GetAntiKt(pythia_particle_list, pythia_true_antikt_jet_0_list, pythia_true_antikt_jet_constituents_0_list, antikt_00);
    GetAntiKt(pythia_particle_list, pythia_true_antikt_jet_1_list, pythia_true_antikt_jet_constituents_1_list, antikt_01);
    GetAntiKt(pythia_particle_list, pythia_true_antikt_jet_2_list, pythia_true_antikt_jet_constituents_2_list, antikt_02);
    GetAntiKt(pythia_particle_list, pythia_true_antikt_jet_3_list, pythia_true_antikt_jet_constituents_3_list, antikt_03);

    DoEmbedding(pythia_true_antikt_jet_constituents_0_list, hijing_matched_jet_subtracted_0_list, embedded_particle_antikt_0_list);
    DoEmbedding(pythia_true_antikt_jet_constituents_1_list, hijing_matched_jet_subtracted_1_list, embedded_particle_antikt_1_list);
    DoEmbedding(pythia_true_antikt_jet_constituents_2_list, hijing_matched_jet_subtracted_2_list, embedded_particle_antikt_2_list);
    DoEmbedding(pythia_true_antikt_jet_constituents_3_list, hijing_matched_jet_subtracted_3_list, embedded_particle_antikt_3_list);

    GetAntiKt(embedded_particle_antikt_0_list,  embedded_good_antikt_jet_0_list, embedded_good_antikt_jet_constituents_0_list, antikt_00);
    GetAntiKt(embedded_particle_antikt_1_list,  embedded_good_antikt_jet_1_list, embedded_good_antikt_jet_constituents_1_list, antikt_01);
    GetAntiKt(embedded_particle_antikt_2_list,  embedded_good_antikt_jet_2_list, embedded_good_antikt_jet_constituents_2_list, antikt_02);
    GetAntiKt(embedded_particle_antikt_3_list,  embedded_good_antikt_jet_3_list, embedded_good_antikt_jet_constituents_3_list, antikt_03);

    GetMatchedJets(pythia_true_antikt_jet_0_list, embedded_good_antikt_jet_0_list, embedded_good_antikt_jet_constituents_0_list,
                   embedded_matched_antikt_jet_0_list, embedded_matched_antikt_jet_constituents_0_list,
                   hDistance_pythia_Antikt_0, hPtTrueJet_pythia_Antikt_0, hPtTrueVsMatched_pythia_Antikt_0,
                   hPtMatchedJet_pythia_Antikt_0, hForMatchingEfficiency_pythia_Antikt_0,
                   minDeltaRAntikt0);
    GetMatchedJets(pythia_true_antikt_jet_1_list, embedded_good_antikt_jet_1_list, embedded_good_antikt_jet_constituents_1_list,
                   embedded_matched_antikt_jet_1_list, embedded_matched_antikt_jet_constituents_1_list,
                   hDistance_pythia_Antikt_1, hPtTrueJet_pythia_Antikt_1, hPtTrueVsMatched_pythia_Antikt_1,
                   hPtMatchedJet_pythia_Antikt_1, hForMatchingEfficiency_pythia_Antikt_1,
                   minDeltaRAntikt1);
    GetMatchedJets(pythia_true_antikt_jet_2_list, embedded_good_antikt_jet_2_list, embedded_good_antikt_jet_constituents_2_list,
                   embedded_matched_antikt_jet_2_list, embedded_matched_antikt_jet_constituents_2_list,
                   hDistance_pythia_Antikt_2, hPtTrueJet_pythia_Antikt_2, hPtTrueVsMatched_pythia_Antikt_2,
                   hPtMatchedJet_pythia_Antikt_2, hForMatchingEfficiency_pythia_Antikt_2,
                   minDeltaRAntikt2);
    GetMatchedJets(pythia_true_antikt_jet_3_list, embedded_good_antikt_jet_3_list, embedded_good_antikt_jet_constituents_3_list,
                   embedded_matched_antikt_jet_3_list, embedded_matched_antikt_jet_constituents_3_list,
                   hDistance_pythia_Antikt_3, hPtTrueJet_pythia_Antikt_3, hPtTrueVsMatched_pythia_Antikt_3,
                   hPtMatchedJet_pythia_Antikt_3, hForMatchingEfficiency_pythia_Antikt_3,
                   minDeltaRAntikt3);

    FillHistograms(embedded_good_antikt_jet_0_list, embedded_matched_antikt_jet_0_list,
                   hPtGoodJet_discriminant_pythia_Antikt_0, hPtMatchedJet_discriminant_pythia_Antikt_0,
                   hPtGoodJet_centrality_pythia_Antikt_0, hPtMatchedJet_centrality_pythia_Antikt_0);
    FillHistograms(embedded_good_antikt_jet_1_list, embedded_matched_antikt_jet_1_list,
                   hPtGoodJet_discriminant_pythia_Antikt_1, hPtMatchedJet_discriminant_pythia_Antikt_1,
                   hPtGoodJet_centrality_pythia_Antikt_1, hPtMatchedJet_centrality_pythia_Antikt_1);
    FillHistograms(embedded_good_antikt_jet_2_list, embedded_matched_antikt_jet_2_list,
                   hPtGoodJet_discriminant_pythia_Antikt_2, hPtMatchedJet_discriminant_pythia_Antikt_2,
                   hPtGoodJet_centrality_pythia_Antikt_2, hPtMatchedJet_centrality_pythia_Antikt_2);
    FillHistograms(embedded_good_antikt_jet_3_list, embedded_matched_antikt_jet_3_list,
                   hPtGoodJet_discriminant_pythia_Antikt_3, hPtMatchedJet_discriminant_pythia_Antikt_3,
                   hPtGoodJet_centrality_pythia_Antikt_3, hPtMatchedJet_centrality_pythia_Antikt_3);

    //*****************************************************************************************************************************
    //Gaussian Filter
    //*****************************************************************************************************************************
    GetGaussianFilter(pythia_particle_list, pythia_true_filter_jet_0_list, pythia_true_filter_jet_constituents_0_list, filter_00, 0.106);
    GetGaussianFilter(pythia_particle_list, pythia_true_filter_jet_1_list, pythia_true_filter_jet_constituents_1_list, filter_01, 0.141);
    GetGaussianFilter(pythia_particle_list, pythia_true_filter_jet_2_list, pythia_true_filter_jet_constituents_2_list, filter_02, 0.177);
    GetGaussianFilter(pythia_particle_list, pythia_true_filter_jet_3_list, pythia_true_filter_jet_constituents_3_list, filter_03, 0.212);

    DoEmbedding(pythia_true_filter_jet_constituents_0_list, hijing_matched_jet_subtracted_0_list, embedded_particle_filter_0_list);
    DoEmbedding(pythia_true_filter_jet_constituents_1_list, hijing_matched_jet_subtracted_1_list, embedded_particle_filter_1_list);
    DoEmbedding(pythia_true_filter_jet_constituents_2_list, hijing_matched_jet_subtracted_2_list, embedded_particle_filter_2_list);
    DoEmbedding(pythia_true_filter_jet_constituents_3_list, hijing_matched_jet_subtracted_3_list, embedded_particle_filter_3_list);

    GetGaussianFilter(embedded_particle_filter_0_list,  embedded_good_filter_jet_0_list,
                      embedded_good_filter_jet_constituents_0_list, filter_00, 0.106);
    GetGaussianFilter(embedded_particle_filter_1_list,  embedded_good_filter_jet_1_list,
                      embedded_good_filter_jet_constituents_1_list, filter_01, 0.141);
    GetGaussianFilter(embedded_particle_filter_2_list,  embedded_good_filter_jet_2_list,
                      embedded_good_filter_jet_constituents_2_list, filter_02, 0.177);
    GetGaussianFilter(embedded_particle_filter_3_list,  embedded_good_filter_jet_3_list,
                      embedded_good_filter_jet_constituents_3_list, filter_03, 0.212);

    GetMatchedJets(pythia_true_filter_jet_0_list, embedded_good_filter_jet_0_list, embedded_good_filter_jet_constituents_0_list,
                   embedded_matched_filter_jet_0_list, embedded_matched_filter_jet_constituents_0_list,
                   hDistance_pythia_Filter_0, hPtTrueJet_pythia_Filter_0, hPtTrueVsMatched_pythia_Filter_0,
                   hPtMatchedJet_pythia_Filter_0, hForMatchingEfficiency_pythia_Filter_0,
                   minDeltaRFilter0);
    GetMatchedJets(pythia_true_filter_jet_1_list, embedded_good_filter_jet_1_list, embedded_good_filter_jet_constituents_1_list,
                   embedded_matched_filter_jet_1_list, embedded_matched_filter_jet_constituents_1_list,
                   hDistance_pythia_Filter_1, hPtTrueJet_pythia_Filter_1, hPtTrueVsMatched_pythia_Filter_1,
                   hPtMatchedJet_pythia_Filter_1, hForMatchingEfficiency_pythia_Filter_1,
                   minDeltaRFilter1);
    GetMatchedJets(pythia_true_filter_jet_2_list, embedded_good_filter_jet_2_list, embedded_good_filter_jet_constituents_2_list,
                   embedded_matched_filter_jet_2_list, embedded_matched_filter_jet_constituents_2_list,
                   hDistance_pythia_Filter_2, hPtTrueJet_pythia_Filter_2, hPtTrueVsMatched_pythia_Filter_2,
                   hPtMatchedJet_pythia_Filter_2, hForMatchingEfficiency_pythia_Filter_2,
                   minDeltaRFilter2);
    GetMatchedJets(pythia_true_filter_jet_3_list, embedded_good_filter_jet_3_list, embedded_good_filter_jet_constituents_3_list,
                   embedded_matched_filter_jet_3_list, embedded_matched_filter_jet_constituents_3_list,
                   hDistance_pythia_Filter_3, hPtTrueJet_pythia_Filter_3, hPtTrueVsMatched_pythia_Filter_3,
                   hPtMatchedJet_pythia_Filter_3, hForMatchingEfficiency_pythia_Filter_3,
                   minDeltaRFilter3);

    FillHistograms(embedded_good_filter_jet_0_list, embedded_matched_filter_jet_0_list,
                   hPtGoodJet_discriminant_pythia_Filter_0, hPtMatchedJet_discriminant_pythia_Filter_0,
                   hPtGoodJet_centrality_pythia_Filter_0, hPtMatchedJet_centrality_pythia_Filter_0);
    FillHistograms(embedded_good_filter_jet_1_list, embedded_matched_filter_jet_1_list,
                   hPtGoodJet_discriminant_pythia_Filter_1, hPtMatchedJet_discriminant_pythia_Filter_1,
                   hPtGoodJet_centrality_pythia_Filter_1, hPtMatchedJet_centrality_pythia_Filter_1);
    FillHistograms(embedded_good_filter_jet_2_list, embedded_matched_filter_jet_2_list,
                   hPtGoodJet_discriminant_pythia_Filter_2, hPtMatchedJet_discriminant_pythia_Filter_2,
                   hPtGoodJet_centrality_pythia_Filter_2, hPtMatchedJet_centrality_pythia_Filter_2);
    FillHistograms(embedded_good_filter_jet_3_list, embedded_matched_filter_jet_3_list,
                   hPtGoodJet_discriminant_pythia_Filter_3, hPtMatchedJet_discriminant_pythia_Filter_3,
                   hPtGoodJet_centrality_pythia_Filter_3, hPtMatchedJet_centrality_pythia_Filter_3);

    //*****************************************************************************************************************************
    //Fill trees
    FillTrees(fillTrees);

    // any other return code might lead to aborting the event or analysis
    return 0;
}

int JetSimWithoutDetector::EndRun(const int runNumber)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }


    if (verbo)
        {

        }
    return 0;
}

int JetSimWithoutDetector::End(PHCompositeNode *topNode)
{
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nRunEvents << endl;

    if(nRunEvents > 0)
        {
            cout << "This is the end. If everthing was ok, you will see this message: " << endl;
            cout << "123ALLDONE321" << endl;
        }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    hEvents->SetBinContent(1, nRunEvents);

    outfile->Write(outfname.c_str());
    outfile->Close();

    delete phpythia;
    delete phhijing;
    delete tpythia6;
    delete accept;

    delete antikt_00;
    delete antikt_01;
    delete antikt_02;
    delete antikt_03;

    delete background_model_perp;
    delete background_model_time;
    delete background_model_z;

    delete filter_00;
    delete filter_01;
    delete filter_02;
    delete filter_03;

    if (verbo)
        {
            cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
        }
    return 0;
}


void JetSimWithoutDetector::GetsHijing()
{
    unsigned int jetId = 0;
    for ( HepMC::GenEvent::particle_iterator p = shijing->particles_begin(); p != shijing->particles_end(); ++p )
        {
            int id = (*p)->pdg_id();

            //No neutrino
            if (id == PY_NU_E || id == -PY_NU_E || id == PY_NU_MU || id == -PY_NU_MU || id == PY_NU_TAU || id == -PY_NU_TAU)
                {
                    continue;
                }

            //No muons
            if (id == PY_MU || id == -PY_MU || id == PY_TAU || id == -PY_TAU)
                {
                    continue;
                }

            const HepMC::FourVector& mom_vector = (*p)->momentum();

            // ignore all beam remnants
            if ( sqrt(pow(mom_vector.px(), 2) + pow(mom_vector.py(), 2)) == 0.0)
                {
                    continue;
                }

            float rho = sqrt(pow(mom_vector.px(), 2) + pow(mom_vector.py(), 2) + pow(mom_vector.pz(), 2));
            float eta = 0.5 * log((rho + mom_vector.pz()) / (rho - mom_vector.pz()));

            float pT = sqrt(pow(mom_vector.px(), 2) + pow(mom_vector.py(), 2));
            float phi = phiReduce(atan2(mom_vector.py(), mom_vector.px()));

            int charge = tpythia6->Pychge(id);
            if (charge == 3)
                {
                    charge = 1;
                }
            if (charge == -3)
                {
                    charge = -1;
                }

            //cut on current PHENIX acceptance
            TLorentzVector *kin = new TLorentzVector(mom_vector.px(), mom_vector.py(), mom_vector.pz(), (*p)->momentum().e());
            float zvertex = 0.0;

            //central acceptance (true = 2pi, false = PHENIX)
            accept->SetFullTwoPi(false);

            bool inPhi, inEta, inBigEta, deadArea;
            int ok = accept->acceptParticle(kin, charge, zvertex, inPhi, inEta, inBigEta, deadArea);
            delete kin;

            bool good = ((ok == 1) || (ok == 2));

            int arm = 1;
            if (phi > 1.57)
                {
                    arm = 0;
                }

            bool acceptanceCut = meetPhenixAcceptanceForTrueJet(eta, phi);

            //Get Fragmentation Jet from sHijing
            if ((*p)->status() == 103 )
                {
                    jets tempJet;
                    tempJet.pT = pT;
                    tempJet.eta = eta;
                    tempJet.phi = phi;
                    tempJet.totalNc = 0.0;
                    tempJet.cf = 0.0;
                    tempJet.nf = 0.0;
                    tempJet.disc = 0.0;
                    tempJet.arm = arm;
                    tempJet.centralityBin = centralityBin;
                    tempJet.id = jetId;

                    if(id == 1500000)
                        {
                            if (acceptanceCut && (pT > minPt))
                                {
                                    hijing_true_antikt_jet_0_list.push_back(tempJet);
                                }
                            continue;
                        }
                    if(id == 2000000)
                        {
                            if (acceptanceCut && (pT > minPt))
                                {
                                    hijing_true_antikt_jet_1_list.push_back(tempJet);
                                }
                            continue;
                        }
                    if(id == 2500000)
                        {
                            if (acceptanceCut && (pT > minPt))
                                {
                                    hijing_true_antikt_jet_2_list.push_back(tempJet);
                                }
                            continue;
                        }
                    if(id == 3000000)
                        {
                            if (acceptanceCut && (pT > minPt))
                                {
                                    hijing_true_antikt_jet_3_list.push_back(tempJet);
                                }
                            continue;
                        }
                }
            jetId++;

            // We want final state particles. end_vertex()- pointer to the decay vertex
            if((*p)->end_vertex() || (*p)->status() != 1 )
                {
                    continue;
                }

            //For centrality
            nHijingMultiplicity++;

            //In PHENIX acceptance
            if (good)
                {
                    particles tempParticle;
                    tempParticle.px = mom_vector.px();
                    tempParticle.py = mom_vector.py();
                    tempParticle.pz = mom_vector.pz();
                    tempParticle.energy = (*p)->momentum().e();
                    tempParticle.pT = pT;
                    tempParticle.eta = eta;
                    tempParticle.phi = phi;
                    tempParticle.mom = sqrt((mom_vector.px() * mom_vector.px()) + (mom_vector.py() * mom_vector.py()) + (mom_vector.pz() * mom_vector.pz()));
                    tempParticle.arm = arm;
                    tempParticle.charge = charge;

                    hijing_particle_list.push_back(tempParticle);
                }
        }
}


void JetSimWithoutDetector::GetPythia()
{
    int npart = phpythia->size();
    for (int ipart = 0; ipart < npart; ipart++)
        {
            TMCParticle *part = phpythia->getParticle(ipart);

            //Look for final state particles in Phenix acceptance
            if (part->GetKS() != 1)
                {
                    continue;
                }

            //KF- falvor
            int kf = part->GetKF();

            //No neutrino
            if (kf == PY_NU_E || kf == -PY_NU_E || kf == PY_NU_MU || kf == -PY_NU_MU || kf == PY_NU_TAU || kf == -PY_NU_TAU)
                {
                    continue;
                }

            //No muons
            if (kf == PY_MU || kf == -PY_MU || kf == PY_TAU || kf == -PY_TAU)
                {
                    continue;
                }

            /**
	       Pythia throws particles with no transverse momentum- and we get this warning:
	       Warning in <TVector3::PseudoRapidity>: transvers momentum = 0! return +/- 10e10
	       -> Throw particles with no pT
            **/
            // ignore all beam remnants
            if ( sqrt(pow(part->GetPx(), 2) + pow(part->GetPy(), 2)) == 0.0)
                {
                    continue;
                }

            float rho = sqrt(pow(part->GetPx(), 2) + pow(part->GetPy(), 2) + pow(part->GetPz(), 2));
            float eta = 0.5 * log((rho + part->GetPz()) / (rho - part->GetPz()));

            float pT = sqrt(pow(part->GetPx(), 2) + pow(part->GetPy(), 2));
            float phi = phiReduce(atan2(part->GetPy(), part->GetPx()));

            //cut on current PHENIX acceptance
            TLorentzVector *kin = new TLorentzVector(part->GetPx(), part->GetPy(), part->GetPz(), part->GetEnergy());
            float zvertex = 0.0;

            int arm = 1;
            if (phi > 1.57)
                {
                    arm = 0;
                }

            int charge = tpythia6->Pychge(kf);
            if (charge == 3)
                {
                    charge = 1;
                }
            if (charge == -3)
                {
                    charge = -1;
                }

            //central acceptance (true = 2pi, false = PHENIX)
            accept->SetFullTwoPi(false);

            bool inPhi, inEta, inBigEta, deadArea;
            int ok = accept->acceptParticle(kin, charge, zvertex, inPhi, inEta, inBigEta, deadArea);
            delete kin;

            bool good = ((ok == 1) || (ok == 2));

            //In PHENIX acceptance
            if (good)
                {
                    particles temp;
                    temp.px = part->GetPx();
                    temp.py = part->GetPy();
                    temp.pz = part->GetPz();
                    temp.energy = part->GetEnergy();
                    temp.pT = pT;
                    temp.eta = eta;
                    temp.phi = phi;
                    temp.mom = sqrt((part->GetPx() * part->GetPx()) + (part->GetPy() * part->GetPy()) + (part->GetPz() * part->GetPz()));
                    temp.arm = arm;
                    temp.charge = charge;

                    pythia_particle_list.push_back(temp);
                }
        }
}

void JetSimWithoutDetector::GetAntiKt(std::vector<particles> particle_list,
                                      std::vector<jets>& antikt_jet_list,
                                      std::vector<constituentParticles>& constituent_list,
                                      fastjet::JetDefinition *antikt,
				      bool trueJet)
{
    //Generate input for jet reconstuction
    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    //Sort the particle list- charged particles first and then neurtal
    std::sort(particle_list.begin(), particle_list.end(), sortParticle());

    unsigned int indexTotal = 0;
    unsigned int indexCharged = 0;
    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            indexTotal++;
            if (particle_list[h].charge != 0)
                {
                    indexCharged++;
                }
        }

    float px[indexTotal];
    float py[indexTotal];
    float pz[indexTotal];
    float energy[indexTotal];
    float particlePt[indexTotal];
    float particleEta[indexTotal];
    float particlePhi[indexTotal];
    float mom[indexTotal];
    int particleArm[indexTotal];
    int charge[indexTotal];

    for (unsigned int cp = 0; cp < indexTotal; cp++)
        {
            px[cp] = 0.0;
            py[cp] = 0.0;
            pz[cp] = 0.0;
            energy[cp] = 0.0;
            particlePt[cp] = 0.0;
            particleEta[cp] = 0.0;
            particlePhi[cp] = 0.0;
            mom[cp] = 0.0;
            particleArm[cp] = 9999;
            charge[cp] = 9999;
        }

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
                    px[index]            = particle_list[h].px;
                    py[index]            = particle_list[h].py;
                    pz[index]            = particle_list[h].pz;
                    energy[index]        = particle_list[h].energy;
                    particlePt[index]    = particle_list[h].pT;
                    particleEta[index]   = particle_list[h].eta;
                    particlePhi[index]   = particle_list[h].phi;
                    mom[index]           = particle_list[h].mom;
                    particleArm[index]   = particle_list[h].arm;
                    charge[index]        = particle_list[h].charge;

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
                    px[index]            = particle_list[h].px;
                    py[index]            = particle_list[h].py;
                    pz[index]            = particle_list[h].pz;
                    energy[index]        = particle_list[h].energy;
                    particlePt[index]    = particle_list[h].pT;
                    particleEta[index]   = particle_list[h].eta;
                    particlePhi[index]   = particle_list[h].phi;
                    mom[index]           = particle_list[h].mom;
                    particleArm[index]   = particle_list[h].arm;
                    charge[index]        = particle_list[h].charge;

                    jetParticles_all.push_back(pseudoNeutral);
                    index++;
                }
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();

    unsigned int jetId = 0;
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float chargedPt = 0.0;
            float neutralPt = 0.0;
            float discriminant = 0.0;

            float jetPt = aFastJet.perp();
            float jetEta = aFastJet.pseudorapidity();
            float jetPhi = phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();
            for (unsigned int iconst = 0; iconst < nconst; iconst++)
                {
                    unsigned int indx = constituents[iconst].user_index();
                    if (indx < indexCharged)//Charged particles
                        {
                            chargedPt += particlePt[indx];
                            float deltaR = dR(particleEta[indx], jetEta, particlePhi[indx], jetPhi);
                            float disc = getDiscriminant(particlePt[indx], deltaR);
                            discriminant += disc;
                        }
                    else //Neutral particles
                        {
                            neutralPt += particlePt[indx];
                            float deltaR = dR(particleEta[indx], jetEta, particlePhi[indx], jetPhi);
                            float disc = getDiscriminant(particlePt[indx], deltaR);
                            discriminant += disc;
                        }
                }

            int arm = 1;
            if (jetPhi > 1.57)
                {
                    arm = 0;
                }

            bool goodAntikt = (jetPt > minPt);

	    if(!trueJet){
		goodAntikt =goodAntiktJet(nconst, jetPt);
	    }

            if (goodAntikt)
                {
                    jets tempJet;
                    tempJet.pT             = jetPt;
                    tempJet.eta            = jetEta;
                    tempJet.phi            = jetPhi;
                    tempJet.totalNc        = (float)nconst;
                    tempJet.cf             = chargedPt / jetPt;
                    tempJet.nf             = neutralPt / jetPt;
                    tempJet.disc           = discriminant;
                    tempJet.arm            = arm;
                    tempJet.centralityBin  = centralityBin;
                    tempJet.id             = jetId;

                    antikt_jet_list.push_back(tempJet);


                    for (unsigned int iconst = 0; iconst < nconst; iconst++)
                        {
                            unsigned int indx = constituents[iconst].user_index();

                            constituentParticles tempConst;
                            tempConst.px      = px[indx];
                            tempConst.py      = py[indx];
                            tempConst.pz      = pz[indx];
                            tempConst.energy  = energy[indx];
                            tempConst.pT      = particlePt[indx];
                            tempConst.eta     = particleEta[indx];
                            tempConst.phi     = particlePhi[indx];
                            tempConst.mom     = mom[indx];
                            tempConst.arm     = particleArm[indx];
                            tempConst.charge  = charge[indx];
                            tempConst.id      = jetId;

                            constituent_list.push_back(tempConst);
                        }
                }
            jetId++;
        }
}

void JetSimWithoutDetector::GetGaussianFilter(std::vector<particles> particle_list,
					      std::vector<jets>& filter_jet_list,
					      std::vector<constituentParticles>& constituent_list,
					      jet::reconstruction_filtering_iir_t *filter,
					      float sigma,
					      bool trueJet)
{
    std::vector<jet::cluster_t> cluster;
    std::vector<jet::track_t> trackAll;

    cluster.clear();
    trackAll.clear();

    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            if (particle_list[h].charge != 0)
                {
                    jet::track_t trackC(NAN, particle_list[h].pT, particle_list[h].eta, particle_list[h].phi);
                    trackAll.push_back(trackC);
                }
            else
                {
                    jet::track_t trackN(NAN, particle_list[h].pT, particle_list[h].eta, particle_list[h].phi);
                    trackAll.push_back(trackN);
                }
        }
    const jet::collision_geometry_t geometry(0.0F, NAN, NAN);

    //All
    jet::event_t test_event_all(0, 0, 1.0F, trackAll, cluster);

    std::vector<jet::jet_t> jet_all;
    jet_all.clear();
    jet_all = filter->reconstruct(test_event_all.track(), geometry);

    unsigned int jetId = 0;
    for (unsigned int n = 0; n < jet_all.size(); n++)
        {
            float totalNc = 0.0;
            float chargedPt = 0.0;
            float neutralPt = 0.0;
            float discriminant = 0.0;

            float gaussPt = jet_all[n].momentum()[1];
            float gaussEta = jet_all[n].momentum()[2];
            float gaussPhi = phiReduce(jet_all[n].momentum()[3]);

            for (unsigned int c = 0; c < particle_list.size(); c++)
                {
                    float hEta = particle_list[c].eta;
                    float hPhi = particle_list[c].phi;
                    float hPt = particle_list[c].pT;

                    float deltaR = dR(hEta, gaussEta, hPhi, gaussPhi);
                    float weight = exp(-deltaR * deltaR / (2.0 * sigma * sigma));
                    float disc = getDiscriminant(hPt, deltaR);
                    if (particle_list[c].charge != 0)
                        {
                            chargedPt +=  hPt * weight;
                            totalNc += weight;
                            discriminant += disc;
                        }
                    else
                        {
                            neutralPt +=  hPt * weight;
                            totalNc += weight;
                            discriminant += disc;
                        }
                }

            int arm = 1;
            if (gaussPhi > 1.57)
                {
                    arm = 0;
                }

	    bool goodFilter = (gaussPt > minPt);

            if(!trueJet){
                goodFilter = goodFilterJet(totalNc, gaussPt);
            }

            if (goodFilter)
                {
                    jets tempJet;
                    tempJet.pT            = gaussPt;
                    tempJet.eta           = gaussEta;
                    tempJet.phi           = gaussPhi;
                    tempJet.totalNc       = totalNc;
                    tempJet.cf            = chargedPt / gaussPt;
                    tempJet.nf            = neutralPt / gaussPt;
                    tempJet.disc          = discriminant;
                    tempJet.arm           = arm;
                    tempJet.centralityBin = centralityBin;
                    tempJet.id            = jetId;

                    filter_jet_list.push_back(tempJet);

                    for (unsigned int c = 0; c < particle_list.size(); c++)
                        {
                            float hEta = particle_list[c].eta;
                            float hPhi = particle_list[c].phi;
                            float hPt = particle_list[c].pT;

                            float deltaR = dR(hEta, gaussEta, hPhi, gaussPhi);
                            float weight = exp(-deltaR * deltaR / (2.0 * sigma * sigma));

                            if (weight > filterWeight)
                                {
                                    constituentParticles tempConst;
                                    tempConst.px     = particle_list[c].px;
                                    tempConst.py     = particle_list[c].py;
                                    tempConst.pz     = particle_list[c].pz;
                                    tempConst.energy = particle_list[c].energy;
                                    tempConst.pT     = particle_list[c].pT;
                                    tempConst.eta    = particle_list[c].eta;
                                    tempConst.phi    = particle_list[c].phi;
                                    tempConst.mom    = particle_list[c].mom;
                                    tempConst.arm    = particle_list[c].arm;
                                    tempConst.charge = particle_list[c].charge;
                                    tempConst.id     = jetId;
                                    constituent_list.push_back(tempConst);
                                }
                        }
                }
            jetId++;
        }
}


void JetSimWithoutDetector::GetMatchedJets(std::vector<jets> true_jet_list,
					   std::vector<jets> good_jet_list,
					   std::vector<constituentParticles> good_constituent_list,
					   std::vector<jets>& matched_jet_list,
					   std::vector<constituentParticles>& matched_constituent_list,
					   TH1D *hDistance, TH1D *hPtTrueJet, TH2D *hPtTrueVsMatched[5],
					   TH1D *hPtMatchedJet, TH1D *hForMatchingEfficiency[5],
					   float minDeltaR)
{
    std::vector<jetPair> jet_pair_list;

    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            float tPt        = true_jet_list[t].pT;
            float tEta       = true_jet_list[t].eta;
            float tPhi       = true_jet_list[t].phi;
            unsigned int tId = true_jet_list[t].id;

            hPtTrueJet->Fill(tPt);

            for(unsigned int g = 0; g < good_jet_list.size(); g++)
                {
                    float gPt        = good_jet_list[g].pT;
                    float gEta       = good_jet_list[g].eta;
                    float gPhi       = good_jet_list[g].phi;
                    unsigned int gId = good_jet_list[g].id;

                    if (true_jet_list[t].arm == good_jet_list[g].arm)
                        {
                            float deltaR = dR(tEta, gEta, tPhi, gPhi);

                            jetPair tempPair;
                            tempPair.id.first = tId;
                            tempPair.id.second = gId;
                            tempPair.pT.first = tPt;
                            tempPair.pT.second = gPt;
                            tempPair.deltaR = deltaR;

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
            jet_pair_list.erase(std::remove_if(jet_pair_list.begin(), jet_pair_list.end(), removePairId(jet_pair_list.front())), jet_pair_list.end());
        }

    //For each unique pair- go through true list and reco list- get true
    for(unsigned int u = 0; u < unique_pair_list.size(); u++)
        {
            unsigned int uniqueIdReco = unique_pair_list[u].id.second;
            float truePt = unique_pair_list[u].pT.first;

            float deltaR = unique_pair_list[u].deltaR;
            hDistance->Fill(deltaR);

            if (deltaR < minDeltaR)
                {
                    for(unsigned int g = 0; g < good_jet_list.size(); g++)
                        {
                            unsigned int gId = good_jet_list[g].id;
                            if(gId == uniqueIdReco)
                                {
                                    float recoPt = good_jet_list[g].pT;

                                    hPtMatchedJet->Fill(recoPt);
                                    matched_jet_list.push_back(good_jet_list[g]);

                                    hPtTrueVsMatched[0]->Fill(recoPt, truePt);
                                    if(centralityBin == 1)
                                        {
                                            hPtTrueVsMatched[1]->Fill(recoPt, truePt);
                                        }
                                    if(centralityBin == 2)
                                        {
                                            hPtTrueVsMatched[2]->Fill(recoPt, truePt);
                                        }
                                    if(centralityBin == 3)
                                        {
                                            hPtTrueVsMatched[3]->Fill(recoPt, truePt);
                                        }
                                    if(centralityBin == 4)
                                        {
                                            hPtTrueVsMatched[4]->Fill(recoPt, truePt);
                                        }

                                    float discriminant = good_jet_list[g].disc;
                                    hForMatchingEfficiency[0]->Fill(truePt);
                                    if (discriminant > discCut1)
                                        {
                                            hForMatchingEfficiency[1]->Fill(truePt);
                                        }
                                    if (discriminant > discCut2)
                                        {
                                            hForMatchingEfficiency[2]->Fill(truePt);
                                        }
                                    if (discriminant > discCut3)
                                        {
                                            hForMatchingEfficiency[3]->Fill(truePt);
                                        }
                                    if (discriminant > discCut4)
                                        {
                                            hForMatchingEfficiency[4]->Fill(truePt);
                                        }
                                }
                        }

                    for(unsigned int gc = 0; gc < good_constituent_list.size(); gc++)
                        {
                            unsigned int gcId = good_constituent_list[gc].id;
                            if(gcId == uniqueIdReco)
                                {
                                    matched_constituent_list.push_back(good_constituent_list[gc]);
                                }
                        }
                }
        }
}

void JetSimWithoutDetector::SubtractMatchedJetConstituents(std::vector<particles> particle_list,
							   std::vector<constituentParticles> matched_constituent_list,
							   std::vector<particles>& hijing_jet_subtracted_list)
{
    for (unsigned int p = 0; p < particle_list.size(); p++)
        {
            particles tempParticle;
            tempParticle.px     = particle_list[p].px;
            tempParticle.py     = particle_list[p].py;
            tempParticle.pz     = particle_list[p].pz;
            tempParticle.energy = particle_list[p].energy;
            tempParticle.pT     = particle_list[p].pT;
            tempParticle.eta    = particle_list[p].eta;
            tempParticle.phi    = particle_list[p].phi;
            tempParticle.mom    = particle_list[p].mom;
            tempParticle.arm    = particle_list[p].arm;
            tempParticle.charge = particle_list[p].charge;

            bool match = false;
            for (unsigned int t = 0; t < matched_constituent_list.size(); t++)
                {
                    int arm = matched_constituent_list[t].arm;
                    int charge = matched_constituent_list[t].charge;
                    float eta = matched_constituent_list[t].eta;
                    float phi = matched_constituent_list[t].phi;
                    float pT = matched_constituent_list[t].pT;
                    float energy = matched_constituent_list[t].energy;

                    bool check = (arm == tempParticle.arm) && (charge == tempParticle.charge) &&
			(fabs(eta - tempParticle.eta) < 0.0001) && (fabs(phi - tempParticle.phi) < 0.0001) &&
			(fabs(pT - tempParticle.pT) < 0.0001) && (fabs(energy - tempParticle.energy) < 0.0001);

                    if (check)
                        {
                            match = true;
                        }
                }

            if (!match)
                {
                    hijing_jet_subtracted_list.push_back(tempParticle);
                }
        }
}

void JetSimWithoutDetector::DoEmbedding(std::vector<constituentParticles> pythia_particle_list,
                                        std::vector<particles> hijing_jet_subtracted_list,
                                        std::vector<particles>& embedded_particle_list)
{
    //Embed particles
    for (unsigned int cp = 0; cp < pythia_particle_list.size(); cp++)
        {
            particles temp;
            temp.px     = pythia_particle_list[cp].px;
            temp.py     = pythia_particle_list[cp].py;
            temp.pz     = pythia_particle_list[cp].pz;
            temp.energy = pythia_particle_list[cp].energy;
            temp.pT     = pythia_particle_list[cp].pT;
            temp.eta    = pythia_particle_list[cp].eta;
            temp.phi    = pythia_particle_list[cp].phi;
            temp.mom    = pythia_particle_list[cp].mom;
            temp.arm    = pythia_particle_list[cp].arm;
            temp.charge = pythia_particle_list[cp].charge;

            embedded_particle_list.push_back(temp);
        }

    for (unsigned int ch = 0; ch < hijing_jet_subtracted_list.size(); ch++)
        {
            particles temp;
            temp.px     = hijing_jet_subtracted_list[ch].px;
            temp.py     = hijing_jet_subtracted_list[ch].py;
            temp.pz     = hijing_jet_subtracted_list[ch].pz;
            temp.energy = hijing_jet_subtracted_list[ch].energy;
            temp.pT     = hijing_jet_subtracted_list[ch].pT;
            temp.eta    = hijing_jet_subtracted_list[ch].eta;
            temp.phi    = hijing_jet_subtracted_list[ch].phi;
            temp.mom    = hijing_jet_subtracted_list[ch].mom;
            temp.arm    = hijing_jet_subtracted_list[ch].arm;
            temp.charge = hijing_jet_subtracted_list[ch].charge;

            embedded_particle_list.push_back(temp);
        }
}

void JetSimWithoutDetector::FillHistograms(std::vector<jets> good_jet_list,
					   std::vector<jets> matched_jet_list,
					   TH1D *hGoodDiscriminant[5],
					   TH1D *hMatchedDiscriminant[5],
					   TH1D *hGoodCentrality[5],
					   TH1D *hMatchedCentrality[5])
{
    for (unsigned int g = 0; g < good_jet_list.size(); g++)
        {
            float pT = good_jet_list[g].pT;
            float discriminant = good_jet_list[g].disc;

            hGoodDiscriminant[0]->Fill(pT);
            if (discriminant > discCut1)
                {
                    hGoodDiscriminant[1]->Fill(pT);
                }
            if (discriminant > discCut2)
                {
                    hGoodDiscriminant[2]->Fill(pT);

                    hGoodCentrality[0]->Fill(pT);
                    if(centralityBin == 1)
                        {
                            hGoodCentrality[1]->Fill(pT);
                        }
                    if(centralityBin == 2)
                        {
                            hGoodCentrality[2]->Fill(pT);
                        }
                    if(centralityBin == 3)
                        {
                            hGoodCentrality[3]->Fill(pT);
                        }
                    if(centralityBin == 4)
                        {
                            hGoodCentrality[4]->Fill(pT);
                        }
                }
            if (discriminant > discCut3)
                {
                    hGoodDiscriminant[3]->Fill(pT);
                }
            if (discriminant > discCut4)
                {
                    hGoodDiscriminant[4]->Fill(pT);
                }
        }

    for (unsigned int g = 0; g < matched_jet_list.size(); g++)
        {
            float pT = matched_jet_list[g].pT;
            float discriminant = matched_jet_list[g].disc;

            hMatchedDiscriminant[0]->Fill(pT);
            if (discriminant > discCut1)
                {
                    hMatchedDiscriminant[1]->Fill(pT);
                }
            if (discriminant > discCut2)
                {
                    hMatchedDiscriminant[2]->Fill(pT);

                    hMatchedCentrality[0]->Fill(pT);
                    if(centralityBin == 1)
                        {
                            hMatchedCentrality[1]->Fill(pT);
                        }
                    if(centralityBin == 2)
                        {
                            hMatchedCentrality[2]->Fill(pT);
                        }
                    if(centralityBin == 3)
                        {
                            hMatchedCentrality[3]->Fill(pT);
                        }
                    if(centralityBin == 4)
                        {
                            hMatchedCentrality[4]->Fill(pT);
                        }
                }
            if (discriminant > discCut3)
                {
                    hMatchedDiscriminant[3]->Fill(pT);
                }
            if (discriminant > discCut4)
                {
                    hMatchedDiscriminant[4]->Fill(pT);
                }
        }
}

void JetSimWithoutDetector::FillTrees(bool fillTrees)
{
    if(fillTrees)
        {
            //*****************************************************************************************************************************
            //sHijing
            //*****************************************************************************************************************************
            finalHijing->Fill();

            hijing_true_antikt_jet_0->Fill();
            hijing_true_antikt_jet_1->Fill();
            hijing_true_antikt_jet_2->Fill();
            hijing_true_antikt_jet_3->Fill();

            hijing_good_antikt_jet_0->Fill();
            hijing_good_antikt_jet_1->Fill();
            hijing_good_antikt_jet_2->Fill();
            hijing_good_antikt_jet_3->Fill();

            hijing_good_antikt_jet_constituents_0->Fill();
            hijing_good_antikt_jet_constituents_1->Fill();
            hijing_good_antikt_jet_constituents_2->Fill();
            hijing_good_antikt_jet_constituents_3->Fill();

            hijing_matched_antikt_jet_0->Fill();
            hijing_matched_antikt_jet_1->Fill();
            hijing_matched_antikt_jet_2->Fill();
            hijing_matched_antikt_jet_3->Fill();

            hijing_matched_antikt_jet_constituents_0->Fill();
            hijing_matched_antikt_jet_constituents_1->Fill();
            hijing_matched_antikt_jet_constituents_2->Fill();
            hijing_matched_antikt_jet_constituents_3->Fill();

            hijing_matched_jet_subtracted_0->Fill();
            hijing_matched_jet_subtracted_1->Fill();
            hijing_matched_jet_subtracted_2->Fill();
            hijing_matched_jet_subtracted_3->Fill();

            //*****************************************************************************************************************************
            //Pythia
            //*****************************************************************************************************************************
            finalPythia->Fill();
            //*****************************************************************************************************************************
            //Anti-kt
            //*****************************************************************************************************************************
            pythia_true_antikt_jet_0->Fill();
            pythia_true_antikt_jet_1->Fill();
            pythia_true_antikt_jet_2->Fill();
            pythia_true_antikt_jet_3->Fill();

            pythia_true_antikt_jet_constituents_0->Fill();
            pythia_true_antikt_jet_constituents_1->Fill();
            pythia_true_antikt_jet_constituents_2->Fill();
            pythia_true_antikt_jet_constituents_3->Fill();

            embedded_particle_antikt_0->Fill();
            embedded_particle_antikt_1->Fill();
            embedded_particle_antikt_2->Fill();
            embedded_particle_antikt_3->Fill();

            embedded_good_antikt_jet_0->Fill();
            embedded_good_antikt_jet_1->Fill();
            embedded_good_antikt_jet_2->Fill();
            embedded_good_antikt_jet_3->Fill();

            embedded_good_antikt_jet_constituents_0->Fill();
            embedded_good_antikt_jet_constituents_1->Fill();
            embedded_good_antikt_jet_constituents_2->Fill();
            embedded_good_antikt_jet_constituents_3->Fill();

            embedded_matched_antikt_jet_0->Fill();
            embedded_matched_antikt_jet_1->Fill();
            embedded_matched_antikt_jet_2->Fill();
            embedded_matched_antikt_jet_3->Fill();

            embedded_matched_antikt_jet_constituents_0->Fill();
            embedded_matched_antikt_jet_constituents_1->Fill();
            embedded_matched_antikt_jet_constituents_2->Fill();
            embedded_matched_antikt_jet_constituents_3->Fill();

            //*****************************************************************************************************************************
            //Gaussian Filter
            //*****************************************************************************************************************************
            pythia_true_filter_jet_0->Fill();
            pythia_true_filter_jet_1->Fill();
            pythia_true_filter_jet_2->Fill();
            pythia_true_filter_jet_3->Fill();

            pythia_true_filter_jet_constituents_0->Fill();
            pythia_true_filter_jet_constituents_1->Fill();
            pythia_true_filter_jet_constituents_2->Fill();
            pythia_true_filter_jet_constituents_3->Fill();

            embedded_particle_filter_0->Fill();
            embedded_particle_filter_1->Fill();
            embedded_particle_filter_2->Fill();
            embedded_particle_filter_3->Fill();

            embedded_good_filter_jet_0->Fill();
            embedded_good_filter_jet_1->Fill();
            embedded_good_filter_jet_2->Fill();
            embedded_good_filter_jet_3->Fill();

            embedded_good_filter_jet_constituents_0->Fill();
            embedded_good_filter_jet_constituents_1->Fill();
            embedded_good_filter_jet_constituents_2->Fill();
            embedded_good_filter_jet_constituents_3->Fill();

            embedded_matched_filter_jet_0->Fill();
            embedded_matched_filter_jet_1->Fill();
            embedded_matched_filter_jet_2->Fill();
            embedded_matched_filter_jet_3->Fill();

            embedded_matched_filter_jet_constituents_0->Fill();
            embedded_matched_filter_jet_constituents_1->Fill();
            embedded_matched_filter_jet_constituents_2->Fill();
            embedded_matched_filter_jet_constituents_3->Fill();
        }
}













