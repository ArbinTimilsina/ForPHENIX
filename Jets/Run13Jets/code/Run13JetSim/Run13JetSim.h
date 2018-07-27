#ifndef __RUN13JETSIM_H__
#define __RUN13JETSIM_H__

#include <Run13JetTriggerPythia.h>
#include <Run13Jet.h>

//Data classes I am using in analysis
#include <Fun4AllServer.h>
#include <recoConsts.h>

const float MM2CM = 0.1;

class Fun4AllHistoManager;
class PHCompositeNode;
class PHPythiaContainer;
class TH1F;
class TFile;
class TTree;

#define NARMSECT 8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72

//pT bins
const int NPTBINS_TRUE = 22;
const float PTBINS_TRUE[NPTBINS_TRUE + 1] =
    {
	5.0000,
	6.0113,
	7.2272,
	8.6890,
	10.4465,
	12.5594,
	15.0998,
	18.1539,
	21.8258,
	26.2404,
	31.5479,
	37.9289,
	45.6005,
	54.8239,
	65.9128,
	79.2447,
	95.2730,
	114.5434,
	137.7114,
	165.5656,
	199.0536,
	239.3150,
	287.7200
    };

const int NPTBINS_RECO = 22;
const float PTBINS_RECO[NPTBINS_RECO + 1] =
    {
	8.0000,
	9.2383,
	10.6682,
	12.3194,
	14.2262,
	16.4282,
	18.9710,
	21.9074,
	25.2982,
	29.2139,
	33.7357,
	38.9574,
	44.9873,
	51.9505,
	59.9915,
	69.2771,
	80.0000,
	92.3826,
	106.6817,
	123.1941,
	142.2624,
	164.2820,
	189.7099
    };


namespace fastjet
{
    class JetDefinition;
}

typedef struct
{
    float pT;
    float eta;
    float phi;
    float nc;
    float cf;

} recoJets;

typedef struct
{
    std::pair<unsigned int, unsigned int> id;
    std::pair<float, float> pT;
    std::pair<float, float> eta;
    std::pair<float, float> phi;
    std::pair<float, float> nc;
    std::pair<float, float> cf;
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


class Run13JetSim: public SubsysReco
{
 public:
    Run13JetSim(float iminPt);
    virtual ~Run13JetSim() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void getPythiaTrueJet(PHCompositeNode *topNode, std::vector<pythiaParticles>& particle_list, std::vector<pythiaJets>& true_jet_list,
                          TH1F *hVertex, TH1F *hPtTrue);

    void getMatchedJet(std::vector<pythiaJets> true_jet_list, std::vector<recoJets> reco_jet_list, std::vector<recoJets>& matched_jet_list,
                       TH1F *hDis, TH1F *hPtTrueM, TH1F *hPtRecoM, TH2F *hResMatrix);

    void fillTrees(bool write);

 protected:
    float minPt;
    bool treesFill;

    TFile* outfile;

    unsigned int nTotalEvents;

    Fun4AllServer* se;
    recoConsts* rc;
    PHPythiaContainer *phpythia;

    Run13JetTriggerPythia *run13JetTriggerPythia;
    Run13Jet *run13Jet;

    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    TH1F *hEvents;
    TH1F *hVertexPythia;
    TH1F *hVertexReco;
    TH1F *hDistance;
    TH1F *hPtTrueJet;
    TH1F *hPtTrueMatchedJet;
    TH1F *hPtRecoJet;
    TH1F *hPtRecoMatchedJet;
    TH2F *hResponseMatrix;

    TH1F *hEMCalTOF;
    TH2F* hSectorHits[NARMSECT];
    TH1F *hClustersEta;
    TH1F *hClustersPhi;

    TH2F *hModifiedQuality_NE;
    TH2F *hModifiedQuality_SE;
    TH2F *hModifiedQuality_NW;
    TH2F *hModifiedQuality_SW;
    TH1F *hTracksZed;
    TH1F *hTracksPhi;
    //*******************************************************************************************************

    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    TTree *tPythiaParticle;
    std::vector<pythiaParticles> pythia_particle;
    TTree *tPythiaTrueJet;
    std::vector<pythiaJets> pythia_true_jet;

    TTree *tRecoTrack;
    std::vector<tracks> charged_particle;
    TTree *tRecoCluster;
    std::vector<clusters> neutral_particle;

    TTree *tRecoJet;
    std::vector<recoJets> reco_jet;
    TTree *tMatchedJet;
    std::vector<recoJets> matched_jet;

    //*******************************************************************************************************
    int getArm(float phi)
    {
        int arm = 1;
        if (phi > 1.57)
            {
                arm = 0;
            }
        return arm;
    }

    float dR(float eta1, float eta2, float phi1, float phi2)
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

    int IsPbGl(int armsect)
    {
        return ((armsect == 4 || armsect == 5) ? 1 : 0);
    }

};

#endif









