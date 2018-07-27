#ifndef __JETPHASESPACE_H__
#define __JETPHASESPACE_H__

#include <JetTriggerPythia.h>
#include <JetAnalyzer.h>

//Data classes I am using in analysis
#include <Fun4AllServer.h>
#include <recoConsts.h>

const float MM2CM = 0.1;
const float minDeltaR = 0.05;
const float minPt = 15.0;

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
        return (jetPair1.id.first == jetPair2.id.first || jetPair1.id.second == jetPair2.id.second);
    }
};

class JetPhaseSpace: public SubsysReco
{
 public:
    JetPhaseSpace();
    virtual ~JetPhaseSpace() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void getPythiaTrueJet(PHCompositeNode *topNode,
                          std::vector<particles>& particle_list,
                          std::vector<jets>& true_jet_list,
                          TH1F *hVertex);

    void getPisaRecoJet(PHCompositeNode *topNode,
                        std::vector<tracks>& all_track_list,
                        std::vector<clusters>& all_cluster_list,
                        std::vector<tracks>& good_track_list,
                        std::vector<clusters>& good_cluster_list,
                        std::vector<particles>& good_particle_list,
                        std::vector<jets>& sim_jet_list,
                        TH1F *hPtReco);

    void getMatchedJet(std::vector<jets> true_jet_list,
                       std::vector<jets> good_jet_list,
                       std::vector<jets>& matched_jet_list,
                       TH1F *hPtTrue,
                       TH1F *hDistance,
                       TH1F *hPtMatched,
                       TH1F *hTotalPhi,
                       TH1F *hPassPhi,
                       TH1F *hTotalEta[5],
                       TH1F *hPassEta[5]);

    void fillTrees(bool what);

 protected:
    bool treesFill;
    TFile* outfile;
    unsigned int nTotalEvents;

    float zVertex;

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
    TH1F *hVertexPythia;
    TH1F *hPtRecoJet;
    TH1F *hDistance;
    TH1F *hPtMatchedJet;

    TH1F *hTotalPhi;
    TH1F *hPassPhi;

    TH1F *hTotalEta[5];
    TH1F *hPassEta[5];

    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************

    TTree *tPythiaParticle;
    std::vector<particles> pythia_particle_list;
    TTree *tPythiaTrueJet;
    std::vector<jets> pythia_true_jet_list;

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
    TTree *tPisaRecoJet;
    std::vector<jets> pisa_reco_jet_list;
    TTree *tPisaRecoMatchedJet;
    std::vector<jets> pisa_reco_matched_jet_list;

    //*******************************************************************************************************

    bool isNominalEta(float zvertex, float eta)
    {
        float thetaPositive = atan2(246, (90 - zvertex));
        if(thetaPositive < -TMath::PiOver2())
            {
                thetaPositive += TMath::TwoPi();
            }
        if(thetaPositive >= 3.0 * TMath::PiOver2())
            {
                thetaPositive -= TMath::TwoPi();
            }

        float thetaNegative = atan2(246, -(90 + zvertex));
        if(thetaNegative < -TMath::PiOver2())
            {
                thetaNegative += TMath::TwoPi();
            }
        if(thetaNegative >= 3.0 * TMath::PiOver2())
            {
                thetaNegative -= TMath::TwoPi();
            }

        float etaEdgeNegative = -log(tan(thetaNegative / 2.0));
        float etaEdgePositive = -log(tan(thetaPositive / 2.0));

        return (eta > etaEdgeNegative && eta < etaEdgePositive);
    }

    bool isNominalPhi(float phi)
    {
        if((phi > DC_PHI_BOTTOM_WEST && phi < DC_PHI_TOP_WEST) || (phi > DC_PHI_TOP_EAST && phi < DC_PHI_BOTTOM_EAST))
            {
                return true;
            }
        else
            {
                return false;
            }
    }
};

#endif
























