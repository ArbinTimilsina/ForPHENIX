#ifndef __FAKEJETSTUDY_H__
#define __FAKEJETSTUDY_H__

#include <JetAnalyzer.h>

#include <Fun4AllServer.h>
#include <recoConsts.h>

class PHCompositeNode;
class PHPythiaContainer;
class TFile;
class TTree;

namespace HepMC
{
    class IO_GenEvent;
    class GenEvent;
};

typedef struct
{
    std::pair<unsigned int, unsigned int> id;
    std::pair<float, float> pT;
    std::pair<float, float> eta;
    std::pair<float, float> phi;
    std::pair<float, float> nc;
    std::pair<float, float> cf;
    std::pair<float, float> nf;
    std::pair<float, float> disc;
    int arm;
    double deltaR;
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
        return (jetPair1.id.second == jetPair2.id.second);
    }
};

const float R             = 0.2;
const float minPt         = MINPT_TRUE;
const float nc            = 3.0;
const float minCf         = 0.2;
const float maxCf         = 0.7;

const float minDeltaR     = 0.2;

class FakeJetStudy: public SubsysReco
{
 public:
    FakeJetStudy();
    virtual ~FakeJetStudy() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);

    void getTrueJet(HepMC::GenEvent *sHijingEvt,
                    std::vector<jets>& true_jet_list,
                    TH1F *hTrue[5]);

    void getRecoJet(PHCompositeNode *topNode,
                    std::vector<tracks>& all_track_list,
                    std::vector<clusters>& all_cluster_list,
                    std::vector<tracks>& good_track_list,
                    std::vector<clusters>& good_cluster_list,
                    std::vector<particles>& good_particle_list,
                    std::vector<jets>& reco_jet_list,
                    std::vector<jets>& reco_jet_for_fake_list,
                    TH1F *hReco[5]);

    void getMatchedJet(std::vector<jets> true_jet_list,
                       std::vector<jets> reco_jet_list,
                       std::vector<jets>& matched_jet_list,
                       std::vector<jets>& not_matched_jet_list,
		       TH1F *hDistance[5],
                       TH2F *hDistanceVsPtTrue[5],
                       TH1F *hTrueMatched[5],
                       TH1F *hRecoMatched[5],
                       TH1F *hRecoNotMatched[5],
                       TH2F *hMatrix[5]);

    void fillTrees(bool what);

 protected:
    TFile* outfile;

    Fun4AllServer* se;
    recoConsts* rc;
    JetAnalyzer *jetAnalyzer;

    bool treesFill;

    unsigned int nTotalEvents;
    unsigned int nTotalEventsCentrality1;
    unsigned int nTotalEventsCentrality2;
    unsigned int nTotalEventsCentrality3;
    unsigned int nTotalEventsCentrality4;

    unsigned int nTrueJetEvents;
    unsigned int nTrueJetEventsCentrality1;
    unsigned int nTrueJetEventsCentrality2;
    unsigned int nTrueJetEventsCentrality3;
    unsigned int nTrueJetEventsCentrality4;

    unsigned int nNoJetEvents;
    unsigned int nNoJetEventsCentrality1;
    unsigned int nNoJetEventsCentrality2;
    unsigned int nNoJetEventsCentrality3;
    unsigned int nNoJetEventsCentrality4;

    int centralityBin;

    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    TH1F *hEvents;

    TH1F *hImpactParameter;
    TH1F *hEventPlane;

    TH1F *hTracksZed[5];
    TH1F *hTracksPhi[5];

    TH1F *hClustersEta[5];
    TH1F *hClustersPhi[5];
    TH1F *hClustersTower[5];

    TH1F *hDistance[5];
    TH2F *hDistanceVsPtTrue[5];
    TH1F *hPtTrue[5];
    TH1F *hPtTrueMatched[5];
    TH1F *hPtReco[5];
    TH1F *hPtRecoMatched[5];
    TH1F *hPtRecoNotMatched[5];
    TH2F *hResponseMatrix[5];

    TH1F *hPtFake[5];
    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    TTree *tTrueJets;
    std::vector<jets> true_jets;

    TTree *tRecoAllTracks;
    std::vector<tracks> reco_all_tracks;
    TTree *tRecoAllClusters;
    std::vector<clusters> reco_all_clusters;
    TTree *tRecoGoodTracks;
    std::vector<tracks> reco_good_tracks;
    TTree *tRecoGoodClusters;
    std::vector<clusters> reco_good_clusters;
    TTree *tRecoGoodParticles;
    std::vector<particles> reco_good_particles;
    TTree *tRecoJets;
    std::vector<jets> reco_jets;
    std::vector<jets> reco_jets_for_fake;

    TTree *tMatchedJets;
    std::vector<jets> matched_jets;
    TTree *tNotMatchedJets;
    std::vector<jets> not_matched_jets;


    //For fake jet
    TTree *tShuffledParticles;
    std::vector<particles> shuffled_particles;
    TTree *tFakeJets;
    std::vector<jets> fake_jets;
    TTree *tFakeConstituents;
    std::vector<particles> fake_jet_constituents;

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

    inline bool meetPhenixAcceptance(float eta, float phi)
    {
	//This is for sHIJING true jet
        return (fabs(eta) < 0.4 && ((phi > -0.64 && phi < 1.03) || (phi > 2.11 && phi < 3.78))); //0.05 away from the edge
    }
};
#endif
























