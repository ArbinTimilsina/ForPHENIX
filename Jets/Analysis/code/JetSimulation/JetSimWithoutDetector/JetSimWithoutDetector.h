#ifndef __JETSIMWITHOUTDETECTOR_H__
#define __JETSIMWITHOUTDETECTOR_H__

#include <SubsysReco.h>
#include <getClass.h>
#include <phool.h>

//Data classes I am using in analysis
#include <PHCentralTrack.h>
#include <PHGlobal.h>
#include <EventHeader.h>
#include <PreviousEvent.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <TriggerHelper.h>
#include <TPythia6.h>
#include <TH1F.h>

//For Acceptance
#include "TAcceptParticle.h"

//For Gaussian Filter
#include <jetrec/rec.h>

const float minDeltaRHijing0 = 0.15;
const float minDeltaRHijing1 = 0.2;
const float minDeltaRHijing2 = 0.25;
const float minDeltaRHijing3 = 0.3;

const float minDeltaRAntikt0 = 0.15;
const float minDeltaRAntikt1 = 0.2;
const float minDeltaRAntikt2 = 0.25;
const float minDeltaRAntikt3 = 0.3;

const float minDeltaRFilter0 = 0.15;
const float minDeltaRFilter1 = 0.2;
const float minDeltaRFilter2 = 0.25;
const float minDeltaRFilter3 = 0.3;

const float minPt = 5.0;
const float filterWeight = 0.4;

const float discCut1 = 10.0;
const float discCut2 = 15.0;
const float discCut3 = 20.0;
const float discCut4 = 25.0;

namespace HepMC
{
    class IO_GenEvent;
    class GenEvent;
};

class Fun4AllHistoManager;
class PHCompositeNode;
class PHPythiaContainer;
class TH1I;
class TH1D;
class TH2D;
class TFile;
class TNtuple;
class TTree;

namespace fastjet
{
    class JetDefinition;
}


typedef struct
{
    float px;
    float py;
    float pz;
    float energy;
    float pT;
    float eta;
    float phi;
    float mom;
    int arm;
    int charge;

} particles;

typedef struct
{
    float px;
    float py;
    float pz;
    float energy;
    float pT;
    float eta;
    float phi;
    float mom;
    int arm;
    int charge;
    unsigned int id;

} constituentParticles;

typedef struct
{
    float pT;
    float eta;
    float phi;
    float totalNc;
    float cf;
    float nf;
    float disc;
    int arm;
    int centralityBin;
    unsigned int id;

} jets;

typedef struct
{
    std::pair<unsigned int, unsigned int> id;
    std::pair<float, float> pT;
    float deltaR;
} jetPair;

struct sortPair
{
    bool operator()(const jetPair &pair1, const jetPair &pair2)
    {
        return (pair1.deltaR < pair2.deltaR);
    }
};

struct removePairId
{
    const jetPair jetPair1;
removePairId(const jetPair &jetPair0)
: jetPair1(jetPair0)
    {}

    bool operator() (const jetPair &jetPair2)
    {
        return (jetPair1.id.first == jetPair2.id.first || jetPair1.id.second == jetPair2.id.second);
    }
};

struct sortParticle
{
    bool operator()(const particles &particle1, const particles &particle2)
    {
        return (fabs(particle1.charge) > fabs(particle2.charge));
    }
};

const int NPTBINS = 25;
const float PTBINS[NPTBINS + 1] =
    {
	5.0000,
	5.4824,
	6.0113,
	6.5913,
	7.2272,
	7.9245,
	8.6890,
	9.5273,
	10.4465,
	11.4543,
	12.5594,
	13.7711,
	15.0998,
	16.5566,
	18.1539,
	19.9054,
	21.8258,
	23.9315,
	26.2404,
	28.7720,
	31.5479,
	34.5915,
	37.9289,
	41.5882,
	45.6005,
	50.0000
    };

class JetSimWithoutDetector: public SubsysReco
{
 public:
    JetSimWithoutDetector(std::string outfilename);
    virtual ~JetSimWithoutDetector() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void GetsHijing();
    void GetPythia();

    void GetAntiKt(std::vector<particles> particle_list,
                   std::vector<jets>& antikt_jet_list,
                   std::vector<constituentParticles>& constituent_list,
                   fastjet::JetDefinition *antikt,
		   bool trueJet = false);

    void GetGaussianFilter(std::vector<particles> particle_list,
                           std::vector<jets>& filter_jet_list,
                           std::vector<constituentParticles>& constituent_list,
                           jet::reconstruction_filtering_iir_t *filter,
                           float sigma,
			   bool trueJet = false);

    void GetMatchedJets(std::vector<jets> true_jet_list,
                        std::vector<jets> good_jet_list,
                        std::vector<constituentParticles> good_constituent_list,
                        std::vector<jets>& matched_jet_list,
                        std::vector<constituentParticles>& matched_constituent_list,
                        TH1D *hDistance, TH1D *hPtTrueJet, TH2D *hPtTrueVsMatched[5],
                        TH1D *hPtMatchedJet, TH1D *hForMatchingEfficiency[5],
                        float minDeltaR);

    void SubtractMatchedJetConstituents(std::vector<particles> particle_list,
                                        std::vector<constituentParticles> matched_constituent_list,
                                        std::vector<particles>& hijing_jet_subtracted_list);

    void DoEmbedding(std::vector<constituentParticles> pythia_particle_list,
                     std::vector<particles> hijing_jet_subtracted_list,
                     std::vector<particles>& embedded_particle_list);

    void FillHistograms(std::vector<jets> good_jet_list,
                        std::vector<jets> matched_jet_list,
                        TH1D *hGoodDiscriminant[5],
                        TH1D *hMatchedDiscriminant[5],
                        TH1D *hGoodCentrality[5],
                        TH1D *hMatchedCentrality[5]);

    void FillTrees(bool fillTrees);

 protected:
    const int verbo;
    const char* hepInput;
    const std::string outfname;
    const bool fillTrees;

    //FastJet: Jet definition
    fastjet::JetDefinition *antikt_00;
    fastjet::JetDefinition *antikt_01;
    fastjet::JetDefinition *antikt_02;
    fastjet::JetDefinition *antikt_03;

    //Gaussian Filter
    jet::factorized_background_model_t *background_model_perp;
    jet::factorized_background_model_t *background_model_time;
    jet::factorized_background_model_t *background_model_z;

    jet::reconstruction_filtering_iir_t *filter_00;
    jet::reconstruction_filtering_iir_t *filter_01;
    jet::reconstruction_filtering_iir_t *filter_02;
    jet::reconstruction_filtering_iir_t *filter_03;

    TFile* outfile;

    PHPythiaContainer *phpythia;
    PHPythiaContainer *phhijing;
    HepMC::GenEvent *shijing;

    TPythia6 *tpythia6;
    TAcceptParticle *accept;

    unsigned int runNumber;
    unsigned int nRunEvents;
    unsigned int nHijingMultiplicity;

    int centralityBin;
    float impactParameter;
    float eventPlane;

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    //General
    //*****************************************************************************************************************************
    TH1F *hEvents;

    TH1D *hForCentrality;
    TH1D *hImpactParameter;
    TH1D *hEventPlane;

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //sHijing
    //*****************************************************************************************************************************
    TH1D *hDistance_hijing_Antikt_0;
    TH1D *hDistance_hijing_Antikt_1;
    TH1D *hDistance_hijing_Antikt_2;
    TH1D *hDistance_hijing_Antikt_3;

    TH1D *hPtTrueJet_hijing_Antikt_0;
    TH1D *hPtTrueJet_hijing_Antikt_1;
    TH1D *hPtTrueJet_hijing_Antikt_2;
    TH1D *hPtTrueJet_hijing_Antikt_3;

    TH1D *hPtMatchedJet_hijing_Antikt_0;
    TH1D *hPtMatchedJet_hijing_Antikt_1;
    TH1D *hPtMatchedJet_hijing_Antikt_2;
    TH1D *hPtMatchedJet_hijing_Antikt_3;

    TH2D *hPtTrueVsMatched_hijing_Antikt_0[5];
    TH2D *hPtTrueVsMatched_hijing_Antikt_1[5];
    TH2D *hPtTrueVsMatched_hijing_Antikt_2[5];
    TH2D *hPtTrueVsMatched_hijing_Antikt_3[5];

    TH1D *hForMatchingEfficiency_hijing_Antikt_0[5];
    TH1D *hForMatchingEfficiency_hijing_Antikt_1[5];
    TH1D *hForMatchingEfficiency_hijing_Antikt_2[5];
    TH1D *hForMatchingEfficiency_hijing_Antikt_3[5];

    TH1D *hPtGoodJet_discriminant_hijing_Antikt_0[5];
    TH1D *hPtGoodJet_discriminant_hijing_Antikt_1[5];
    TH1D *hPtGoodJet_discriminant_hijing_Antikt_2[5];
    TH1D *hPtGoodJet_discriminant_hijing_Antikt_3[5];

    TH1D *hPtMatchedJet_discriminant_hijing_Antikt_0[5];
    TH1D *hPtMatchedJet_discriminant_hijing_Antikt_1[5];
    TH1D *hPtMatchedJet_discriminant_hijing_Antikt_2[5];
    TH1D *hPtMatchedJet_discriminant_hijing_Antikt_3[5];

    TH1D *hPtGoodJet_centrality_hijing_Antikt_0[5];
    TH1D *hPtGoodJet_centrality_hijing_Antikt_1[5];
    TH1D *hPtGoodJet_centrality_hijing_Antikt_2[5];
    TH1D *hPtGoodJet_centrality_hijing_Antikt_3[5];

    TH1D *hPtMatchedJet_centrality_hijing_Antikt_0[5];
    TH1D *hPtMatchedJet_centrality_hijing_Antikt_1[5];
    TH1D *hPtMatchedJet_centrality_hijing_Antikt_2[5];
    TH1D *hPtMatchedJet_centrality_hijing_Antikt_3[5];

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Pythia
    //*****************************************************************************************************************************
    //Anti-kt
    //*****************************************************************************************************************************
    TH1D *hDistance_pythia_Antikt_0;
    TH1D *hDistance_pythia_Antikt_1;
    TH1D *hDistance_pythia_Antikt_2;
    TH1D *hDistance_pythia_Antikt_3;

    TH1D *hPtTrueJet_pythia_Antikt_0;
    TH1D *hPtTrueJet_pythia_Antikt_1;
    TH1D *hPtTrueJet_pythia_Antikt_2;
    TH1D *hPtTrueJet_pythia_Antikt_3;

    TH1D *hPtMatchedJet_pythia_Antikt_0;
    TH1D *hPtMatchedJet_pythia_Antikt_1;
    TH1D *hPtMatchedJet_pythia_Antikt_2;
    TH1D *hPtMatchedJet_pythia_Antikt_3;

    TH2D *hPtTrueVsMatched_pythia_Antikt_0[5];
    TH2D *hPtTrueVsMatched_pythia_Antikt_1[5];
    TH2D *hPtTrueVsMatched_pythia_Antikt_2[5];
    TH2D *hPtTrueVsMatched_pythia_Antikt_3[5];

    TH1D *hForMatchingEfficiency_pythia_Antikt_0[5];
    TH1D *hForMatchingEfficiency_pythia_Antikt_1[5];
    TH1D *hForMatchingEfficiency_pythia_Antikt_2[5];
    TH1D *hForMatchingEfficiency_pythia_Antikt_3[5];

    TH1D *hPtGoodJet_discriminant_pythia_Antikt_0[5];
    TH1D *hPtGoodJet_discriminant_pythia_Antikt_1[5];
    TH1D *hPtGoodJet_discriminant_pythia_Antikt_2[5];
    TH1D *hPtGoodJet_discriminant_pythia_Antikt_3[5];

    TH1D *hPtMatchedJet_discriminant_pythia_Antikt_0[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Antikt_1[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Antikt_2[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Antikt_3[5];

    TH1D *hPtGoodJet_centrality_pythia_Antikt_0[5];
    TH1D *hPtGoodJet_centrality_pythia_Antikt_1[5];
    TH1D *hPtGoodJet_centrality_pythia_Antikt_2[5];
    TH1D *hPtGoodJet_centrality_pythia_Antikt_3[5];

    TH1D *hPtMatchedJet_centrality_pythia_Antikt_0[5];
    TH1D *hPtMatchedJet_centrality_pythia_Antikt_1[5];
    TH1D *hPtMatchedJet_centrality_pythia_Antikt_2[5];
    TH1D *hPtMatchedJet_centrality_pythia_Antikt_3[5];

    //*****************************************************************************************************************************
    //Gaussian Filter
    //*****************************************************************************************************************************
    TH1D *hDistance_pythia_Filter_0;
    TH1D *hDistance_pythia_Filter_1;
    TH1D *hDistance_pythia_Filter_2;
    TH1D *hDistance_pythia_Filter_3;

    TH1D *hPtTrueJet_pythia_Filter_0;
    TH1D *hPtTrueJet_pythia_Filter_1;
    TH1D *hPtTrueJet_pythia_Filter_2;
    TH1D *hPtTrueJet_pythia_Filter_3;

    TH1D *hPtMatchedJet_pythia_Filter_0;
    TH1D *hPtMatchedJet_pythia_Filter_1;
    TH1D *hPtMatchedJet_pythia_Filter_2;
    TH1D *hPtMatchedJet_pythia_Filter_3;

    TH2D *hPtTrueVsMatched_pythia_Filter_0[5];
    TH2D *hPtTrueVsMatched_pythia_Filter_1[5];
    TH2D *hPtTrueVsMatched_pythia_Filter_2[5];
    TH2D *hPtTrueVsMatched_pythia_Filter_3[5];

    TH1D *hForMatchingEfficiency_pythia_Filter_0[5];
    TH1D *hForMatchingEfficiency_pythia_Filter_1[5];
    TH1D *hForMatchingEfficiency_pythia_Filter_2[5];
    TH1D *hForMatchingEfficiency_pythia_Filter_3[5];

    TH1D *hPtGoodJet_discriminant_pythia_Filter_0[5];
    TH1D *hPtGoodJet_discriminant_pythia_Filter_1[5];
    TH1D *hPtGoodJet_discriminant_pythia_Filter_2[5];
    TH1D *hPtGoodJet_discriminant_pythia_Filter_3[5];

    TH1D *hPtMatchedJet_discriminant_pythia_Filter_0[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Filter_1[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Filter_2[5];
    TH1D *hPtMatchedJet_discriminant_pythia_Filter_3[5];

    TH1D *hPtGoodJet_centrality_pythia_Filter_0[5];
    TH1D *hPtGoodJet_centrality_pythia_Filter_1[5];
    TH1D *hPtGoodJet_centrality_pythia_Filter_2[5];
    TH1D *hPtGoodJet_centrality_pythia_Filter_3[5];

    TH1D *hPtMatchedJet_centrality_pythia_Filter_0[5];
    TH1D *hPtMatchedJet_centrality_pythia_Filter_1[5];
    TH1D *hPtMatchedJet_centrality_pythia_Filter_2[5];
    TH1D *hPtMatchedJet_centrality_pythia_Filter_3[5];

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Ntuples
    //*****************************************************************************************************************************
    //sHijing
    //*****************************************************************************************************************************
    TTree *finalHijing;
    std::vector<particles> hijing_particle_list;

    TTree *hijing_true_antikt_jet_0;
    TTree *hijing_true_antikt_jet_1;
    TTree *hijing_true_antikt_jet_2;
    TTree *hijing_true_antikt_jet_3;
    std::vector<jets> hijing_true_antikt_jet_0_list;
    std::vector<jets> hijing_true_antikt_jet_1_list;
    std::vector<jets> hijing_true_antikt_jet_2_list;
    std::vector<jets> hijing_true_antikt_jet_3_list;

    TTree *hijing_good_antikt_jet_0;
    TTree *hijing_good_antikt_jet_1;
    TTree *hijing_good_antikt_jet_2;
    TTree *hijing_good_antikt_jet_3;
    std::vector<jets> hijing_good_antikt_jet_0_list;
    std::vector<jets> hijing_good_antikt_jet_1_list;
    std::vector<jets> hijing_good_antikt_jet_2_list;
    std::vector<jets> hijing_good_antikt_jet_3_list;

    TTree *hijing_good_antikt_jet_constituents_0;
    TTree *hijing_good_antikt_jet_constituents_1;
    TTree *hijing_good_antikt_jet_constituents_2;
    TTree *hijing_good_antikt_jet_constituents_3;
    std::vector<constituentParticles> hijing_good_antikt_jet_constituents_0_list;
    std::vector<constituentParticles> hijing_good_antikt_jet_constituents_1_list;
    std::vector<constituentParticles> hijing_good_antikt_jet_constituents_2_list;
    std::vector<constituentParticles> hijing_good_antikt_jet_constituents_3_list;

    TTree *hijing_matched_antikt_jet_0;
    TTree *hijing_matched_antikt_jet_1;
    TTree *hijing_matched_antikt_jet_2;
    TTree *hijing_matched_antikt_jet_3;
    std::vector<jets> hijing_matched_antikt_jet_0_list;
    std::vector<jets> hijing_matched_antikt_jet_1_list;
    std::vector<jets> hijing_matched_antikt_jet_2_list;
    std::vector<jets> hijing_matched_antikt_jet_3_list;

    TTree *hijing_matched_antikt_jet_constituents_0;
    TTree *hijing_matched_antikt_jet_constituents_1;
    TTree *hijing_matched_antikt_jet_constituents_2;
    TTree *hijing_matched_antikt_jet_constituents_3;
    std::vector<constituentParticles> hijing_matched_antikt_jet_constituents_0_list;
    std::vector<constituentParticles> hijing_matched_antikt_jet_constituents_1_list;
    std::vector<constituentParticles> hijing_matched_antikt_jet_constituents_2_list;
    std::vector<constituentParticles> hijing_matched_antikt_jet_constituents_3_list;

    TTree *hijing_matched_jet_subtracted_0;
    TTree *hijing_matched_jet_subtracted_1;
    TTree *hijing_matched_jet_subtracted_2;
    TTree *hijing_matched_jet_subtracted_3;
    std::vector<particles> hijing_matched_jet_subtracted_0_list;
    std::vector<particles> hijing_matched_jet_subtracted_1_list;
    std::vector<particles> hijing_matched_jet_subtracted_2_list;
    std::vector<particles> hijing_matched_jet_subtracted_3_list;

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Pythia
    TTree *finalPythia;
    std::vector<particles> pythia_particle_list;

    //*****************************************************************************************************************************
    //Anti-kt
    //*****************************************************************************************************************************
    TTree *pythia_true_antikt_jet_0;
    TTree *pythia_true_antikt_jet_1;
    TTree *pythia_true_antikt_jet_2;
    TTree *pythia_true_antikt_jet_3;
    std::vector<jets> pythia_true_antikt_jet_0_list;
    std::vector<jets> pythia_true_antikt_jet_1_list;
    std::vector<jets> pythia_true_antikt_jet_2_list;
    std::vector<jets> pythia_true_antikt_jet_3_list;

    TTree *pythia_true_antikt_jet_constituents_0;
    TTree *pythia_true_antikt_jet_constituents_1;
    TTree *pythia_true_antikt_jet_constituents_2;
    TTree *pythia_true_antikt_jet_constituents_3;
    std::vector<constituentParticles> pythia_true_antikt_jet_constituents_0_list;
    std::vector<constituentParticles> pythia_true_antikt_jet_constituents_1_list;
    std::vector<constituentParticles> pythia_true_antikt_jet_constituents_2_list;
    std::vector<constituentParticles> pythia_true_antikt_jet_constituents_3_list;

    TTree *embedded_particle_antikt_0;
    TTree *embedded_particle_antikt_1;
    TTree *embedded_particle_antikt_2;
    TTree *embedded_particle_antikt_3;
    std::vector<particles> embedded_particle_antikt_0_list;
    std::vector<particles> embedded_particle_antikt_1_list;
    std::vector<particles> embedded_particle_antikt_2_list;
    std::vector<particles> embedded_particle_antikt_3_list;

    TTree *embedded_good_antikt_jet_0;
    TTree *embedded_good_antikt_jet_1;
    TTree *embedded_good_antikt_jet_2;
    TTree *embedded_good_antikt_jet_3;
    std::vector<jets> embedded_good_antikt_jet_0_list;
    std::vector<jets> embedded_good_antikt_jet_1_list;
    std::vector<jets> embedded_good_antikt_jet_2_list;
    std::vector<jets> embedded_good_antikt_jet_3_list;

    TTree *embedded_good_antikt_jet_constituents_0;
    TTree *embedded_good_antikt_jet_constituents_1;
    TTree *embedded_good_antikt_jet_constituents_2;
    TTree *embedded_good_antikt_jet_constituents_3;
    std::vector<constituentParticles> embedded_good_antikt_jet_constituents_0_list;
    std::vector<constituentParticles> embedded_good_antikt_jet_constituents_1_list;
    std::vector<constituentParticles> embedded_good_antikt_jet_constituents_2_list;
    std::vector<constituentParticles> embedded_good_antikt_jet_constituents_3_list;

    TTree *embedded_matched_antikt_jet_0;
    TTree *embedded_matched_antikt_jet_1;
    TTree *embedded_matched_antikt_jet_2;
    TTree *embedded_matched_antikt_jet_3;
    std::vector<jets> embedded_matched_antikt_jet_0_list;
    std::vector<jets> embedded_matched_antikt_jet_1_list;
    std::vector<jets> embedded_matched_antikt_jet_2_list;
    std::vector<jets> embedded_matched_antikt_jet_3_list;

    TTree *embedded_matched_antikt_jet_constituents_0;
    TTree *embedded_matched_antikt_jet_constituents_1;
    TTree *embedded_matched_antikt_jet_constituents_2;
    TTree *embedded_matched_antikt_jet_constituents_3;
    std::vector<constituentParticles> embedded_matched_antikt_jet_constituents_0_list;
    std::vector<constituentParticles> embedded_matched_antikt_jet_constituents_1_list;
    std::vector<constituentParticles> embedded_matched_antikt_jet_constituents_2_list;
    std::vector<constituentParticles> embedded_matched_antikt_jet_constituents_3_list;

    //*****************************************************************************************************************************
    //Gaussian Filter
    //*****************************************************************************************************************************
    TTree *pythia_true_filter_jet_0;
    TTree *pythia_true_filter_jet_1;
    TTree *pythia_true_filter_jet_2;
    TTree *pythia_true_filter_jet_3;
    std::vector<jets> pythia_true_filter_jet_0_list;
    std::vector<jets> pythia_true_filter_jet_1_list;
    std::vector<jets> pythia_true_filter_jet_2_list;
    std::vector<jets> pythia_true_filter_jet_3_list;

    TTree *pythia_true_filter_jet_constituents_0;
    TTree *pythia_true_filter_jet_constituents_1;
    TTree *pythia_true_filter_jet_constituents_2;
    TTree *pythia_true_filter_jet_constituents_3;
    std::vector<constituentParticles> pythia_true_filter_jet_constituents_0_list;
    std::vector<constituentParticles> pythia_true_filter_jet_constituents_1_list;
    std::vector<constituentParticles> pythia_true_filter_jet_constituents_2_list;
    std::vector<constituentParticles> pythia_true_filter_jet_constituents_3_list;

    TTree *embedded_particle_filter_0;
    TTree *embedded_particle_filter_1;
    TTree *embedded_particle_filter_2;
    TTree *embedded_particle_filter_3;
    std::vector<particles> embedded_particle_filter_0_list;
    std::vector<particles> embedded_particle_filter_1_list;
    std::vector<particles> embedded_particle_filter_2_list;
    std::vector<particles> embedded_particle_filter_3_list;

    TTree *embedded_good_filter_jet_0;
    TTree *embedded_good_filter_jet_1;
    TTree *embedded_good_filter_jet_2;
    TTree *embedded_good_filter_jet_3;
    std::vector<jets> embedded_good_filter_jet_0_list;
    std::vector<jets> embedded_good_filter_jet_1_list;
    std::vector<jets> embedded_good_filter_jet_2_list;
    std::vector<jets> embedded_good_filter_jet_3_list;

    TTree *embedded_good_filter_jet_constituents_0;
    TTree *embedded_good_filter_jet_constituents_1;
    TTree *embedded_good_filter_jet_constituents_2;
    TTree *embedded_good_filter_jet_constituents_3;
    std::vector<constituentParticles> embedded_good_filter_jet_constituents_0_list;
    std::vector<constituentParticles> embedded_good_filter_jet_constituents_1_list;
    std::vector<constituentParticles> embedded_good_filter_jet_constituents_2_list;
    std::vector<constituentParticles> embedded_good_filter_jet_constituents_3_list;

    TTree *embedded_matched_filter_jet_0;
    TTree *embedded_matched_filter_jet_1;
    TTree *embedded_matched_filter_jet_2;
    TTree *embedded_matched_filter_jet_3;
    std::vector<jets> embedded_matched_filter_jet_0_list;
    std::vector<jets> embedded_matched_filter_jet_1_list;
    std::vector<jets> embedded_matched_filter_jet_2_list;
    std::vector<jets> embedded_matched_filter_jet_3_list;

    TTree *embedded_matched_filter_jet_constituents_0;
    TTree *embedded_matched_filter_jet_constituents_1;
    TTree *embedded_matched_filter_jet_constituents_2;
    TTree *embedded_matched_filter_jet_constituents_3;
    std::vector<constituentParticles> embedded_matched_filter_jet_constituents_0_list;
    std::vector<constituentParticles> embedded_matched_filter_jet_constituents_1_list;
    std::vector<constituentParticles> embedded_matched_filter_jet_constituents_2_list;
    std::vector<constituentParticles> embedded_matched_filter_jet_constituents_3_list;

    inline float dR(float eta1, float eta2, float phi1, float phi2)
    {
        float deta = eta1 - eta2;
        float dphi = phi1 - phi2;
        if (dphi < -TMath::Pi())
            {
                dphi += 2 * TMath::Pi();
            }
        if (dphi > TMath::Pi())
            {
                dphi -= 2 * TMath::Pi();
            }
        return sqrt(dphi * dphi + deta * deta);
    }

    inline float dPhi(float phi1, float phi2)
    {
        float dphi = phi1 - phi2;
        if (dphi < -TMath::Pi())
            {
                dphi += 2 * TMath::Pi();
            }
        if (dphi > TMath::Pi())
            {
                dphi -= 2 * TMath::Pi();
            }
        return dphi;
    }

    inline float phiReduce(float phi)
    {
        if (phi < -TMath::PiOver2())
            {
                phi += TMath::TwoPi();
            }
        if (phi >= 3.0 * TMath::PiOver2())
            {
                phi -= TMath::TwoPi();
            }
        return phi;
    }

    inline float getDiscriminant(float pT, float dR)
    {
        return (pT * pT * exp(-dR * dR / (2.0 * 0.1 * 0.1)));
    }

    inline bool meetFiducialCut(float eta, float phi)
    {
        return (fabs(eta) < 0.30 && ((phi > -0.54 && phi < 0.93) || (phi > 2.21 && phi < 3.68))); //0.05 from the edge
    }

    inline bool meetPhenixAcceptance(float eta, float phi)
    {
        return (fabs(eta) < 0.35 && ((phi > -0.59 && phi < 0.98) || (phi > 2.16 && phi < 3.73)));
    }

    inline bool meetPhenixAcceptanceForTrueJet(float eta, float phi)
    {
        return (fabs(eta) < 0.4 && ((phi > -0.64 && phi < 1.03) || (phi > 2.11 && phi < 3.78))); //0.05 away from the edge
    }

    inline bool goodAntiktJet(int nconst, float pT)
    {
        return ((nconst >= 3) && (pT > minPt));
    }

    inline bool goodFilterJet(float nconst, float pT)
    {
        return ((nconst >= 2.5) && (pT > minPt));
    }

    inline int getPtBin(float pT)
    {
        int pbin = -1;
        if (pT >= 5.0 && pT < 7.0)
            {
                pbin = 0;
            }
        if (pT >= 7.0 && pT < 10.0)
            {
                pbin = 1;
            }
        if (pT >= 10.0)
            {
                pbin = 2;
            }
        return pbin;
    }

    inline int getCentralityBin(float iP)
    {
        int cbin = -1;
        if(iP > 0.0 && iP <= 5.634)
            {
                cbin = 1;
            }
        if(iP > 5.634 && iP <= 8.085)
            {
                cbin = 2;
            }
        if(iP > 8.085 && iP <= 10.04)
            {
                cbin = 3;
            }
        if(iP > 10.04 && iP <= 12.95)
            {
                cbin = 4;
            }
        return cbin;
    }
};


#endif








