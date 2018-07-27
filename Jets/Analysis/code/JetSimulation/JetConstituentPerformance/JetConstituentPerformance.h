#ifndef __JETCONSTITUENTPERFORMANCE_H__
#define __JETCONSTITUENTPERFORMANCE_H__

#include <JetTriggerPythia.h>
#include <JetAnalyzer.h>

//Data classes I am using in analysis
#include <Fun4AllServer.h>
#include <recoConsts.h>

class Fun4AllHistoManager;
class PHCompositeNode;
class PHPythiaContainer;
class TH1F;
class TFile;
class TTree;

namespace fastjet
{
    class JetDefinition;
}

typedef struct
{
    std::pair<unsigned int, unsigned int> id;
    std::pair<float, float> pT;
    std::pair<float, float> eT;
    double deltaR;
} particlePair;

struct sortPair
{
    bool operator()(const particlePair &pair1, const particlePair &pair2)
    {
        return (pair1.deltaR < pair2.deltaR);
    }
};

struct removePairId
{
    const particlePair particlePair1;
removePairId(const particlePair &particlePair0)
: particlePair1(particlePair0)
    {}

    bool operator() (const particlePair &particlePair2)
    {
        return (particlePair1.id.first == particlePair2.id.first || particlePair1.id.second == particlePair2.id.second);
    }
};

class JetConstituentPerformance: public SubsysReco
{
 public:
    JetConstituentPerformance(const float ir, const float iminPt,
                              const float incPythia, const float iminCfPythia, const float imaxCfPythia,
                              const float iminChargedDeltaR, const float iminNeutralDeltaR,
                              const bool iperfectEMCal = false, const bool iphenixParticle = false);

    virtual ~JetConstituentPerformance() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);
    void SetCuAu     (bool what);

    void getPythiaTrueJet(PHCompositeNode *topNode,
                          std::vector<particles>& particle_list,
                          std::vector<jets>& true_jet_list,
                          std::vector<particles>& true_jet_constituent_list);

    void getEfficiency(PHCompositeNode *topNode,
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
                       bool perfectEMC);

    void getTrackCutsEfficiency(PHCompositeNode *topNode,
                                std::vector<particles> true_jet_constituent_list,
                                TH1F *hTDistance[4],
                                TH1F *hTEfficiency[4]);

    void fillHistograms();
    void fillTrees(bool what);

 protected:
    float R;
    float minPt;
    float ncPythia;
    float minCfPythia;
    float maxCfPythia;

    float minChargedDeltaR;
    float minNeutralDeltaR;

    bool perfectEMCal;
    bool phenixParticle;

    bool treesFill;

    TFile* outfile;

    unsigned int nTotalEvents;
    unsigned int nTotalPyParticles;

    Fun4AllServer* se;
    recoConsts* rc;
    PHPythiaContainer *phpythia;
    JetAnalyzer *jetAnalyzer;
    JetTriggerPythia *jetTriggerPythia;

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    TH1F *hEvents;

    //*******************************************************************************************************
    TH1F *hPtTrueJet;
    TH1F *hPtPion;
    TH1F *hPtPhoton;
    TH1F *hPtKaon;
    TH1F *hPtProton;
    TH1F *hPtElectron;
    TH1F *hPtNeutron;
    TH1F *hPtK0L;

    //*******************************************************************************************************
    TH1F *hPEtTrueParticle[2];
    TH1F *hDistance[2];
    TH1F *hPEtRecoParticle[2];
    TH1F *hPEtForEfficiency[2];
    TH2F *hResponseMatrix[2];

    //*******************************************************************************************************
    TH1F *hTracksDistance[4];
    TH1F *hTracksEfficiency[4];

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************

    TTree *tPythiaParticle;
    std::vector<particles> pythia_particle_list;
    TTree *tPythiaTrueJet;
    std::vector<jets> pythia_true_jet_list;
    TTree *tPythiaTrueJetConstituent;
    std::vector<particles> pythia_true_jet_constituent_list;

    //*******************************************************************************************************
    TTree *tPisaRecoAllTrack;
    std::vector<tracks> pisa_reco_all_track_list;
    TTree *tPisaRecoAllCluster;
    std::vector<clusters> pisa_reco_all_cluster_list;

    TTree *tPisaRecoGoodTrack;
    std::vector<tracks> pisa_reco_good_track_list;
    TTree *tPisaRecoGoodCluster;
    std::vector<clusters> pisa_reco_good_cluster_list;
    TTree *tPisaRecoGoodParticle;
    std::vector<particles> pisa_reco_good_particle_list;

    //*******************************************************************************************************
};

#endif





















