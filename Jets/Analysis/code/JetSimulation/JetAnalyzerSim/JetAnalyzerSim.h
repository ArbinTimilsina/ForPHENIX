#ifndef __JETANALYZERSIM_H__
#define __JETANALYZERSIM_H__

#include <JetAnalyzer.h>

//Data classes I am using in analysis
#include <Fun4AllServer.h>
#include <recoConsts.h>

const float MM2CM = 0.1;

class Fun4AllHistoManager;
class PHCompositeNode;
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
    int centralityBin;
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

class JetAnalyzerSim: public SubsysReco
{
 public:
    JetAnalyzerSim(const float ir, const float iminPt,
                   const float inc, const float iminCf, const float imaxCf,
                   const float iminDeltaR);

    virtual ~JetAnalyzerSim() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void getPythiaTrueJet(std::vector<jets> true_jet_list,
                          TH1F *hPtTrue[5],
                          TH1F *hPtTrueEast[5],
                          TH1F *hPtTrueWest[5]);

    void getRecoJet(std::vector<particles> total_particle_list,
                    std::vector<jets>& reco_jet_list,
                    TH1F *hPtReco[5]);

    void getMatchedJet(float fMinDeltaR,
                       std::vector<jets> true_jet_list,
                       std::vector<jets> reco_jet_list,
                       std::vector<jets>& matched_jet_list,
                       TH1F *hDistance[5],
                       TH2F *hDistanceVsPtTrue[5],
                       TH1F *hPtTrueMatched[5],
                       TH1F *hPtRecoMatched[5],
                       TH2F *hResponseMatrix[5],
                       TH1F *hForJESR[5][22]);

    void getVariation(float vNc, float vMinCf, float vMaxCf, int fiducialCut,
                      std::vector<jets> true_jet_list,
                      std::vector<particles> total_particle_list,
                      TH1F *hPtReco[5],
                      TH1F *hPtTrueMatched[5],
                      TH1F *hPtRecoMatched[5],
                      TH2F *hResponseMatrix[5]);

    void fillTrees(bool what);

 protected:
    float R;
    float minPt;
    float nc;
    float minCf;
    float maxCf;

    float minDeltaR;

    bool treesFill;

    int centralityBin;

    TFile* outfile;

    JetAnalyzer *jetAnalyzer;

    unsigned int nTotalEvents;
    unsigned int nTotalEventsCentrality1;
    unsigned int nTotalEventsCentrality2;
    unsigned int nTotalEventsCentrality3;
    unsigned int nTotalEventsCentrality4;

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    TH1F *hEvents;

    //*******************************************************************************************************

    TH1F *hPtTrueJets[5];
    TH1F *hPtTrueJetsEast[5];
    TH1F *hPtTrueJetsWest[5];

    //*******************************************************************************************************

    TH1F *hPtRecoJets[5];

    TH1F *hDistance[5];
    TH2F *hDistanceVsPtTrue[5];
    TH1F *hPtTrueMatchedJets[5];
    TH1F *hPtRecoMatchedJets[5];
    TH2F *hResponseMatrix[5];
    TH1F *hForJESR[5][22];

    //*******************************************************************************************************
    //Variations
    //*******************************************************************************************************
    TH1F *hPtRecoDefault[5];
    TH1F *hPtTrueMatchedDefault[5];
    TH1F *hPtRecoMatchedDefault[5];
    TH2F *hResponseMatrixDefault[5];

    TH1F *hPtRecoEast[5];
    TH1F *hPtTrueMatchedEast[5];
    TH1F *hPtRecoMatchedEast[5];
    TH2F *hResponseMatrixEast[5];

    TH1F *hPtRecoWest[5];
    TH1F *hPtTrueMatchedWest[5];
    TH1F *hPtRecoMatchedWest[5];
    TH2F *hResponseMatrixWest[5];

    TH1F *hPtRecoFidTight[5];
    TH1F *hPtTrueMatchedFidTight[5];
    TH1F *hPtRecoMatchedFidTight[5];
    TH2F *hResponseMatrixFidTight[5];

    TH1F *hPtRecoNc[5];
    TH1F *hPtTrueMatchedNc[5];
    TH1F *hPtRecoMatchedNc[5];
    TH2F *hResponseMatrixNc[5];

    TH1F *hPtRecoCf[5];
    TH1F *hPtTrueMatchedCf[5];
    TH1F *hPtRecoMatchedCf[5];
    TH2F *hResponseMatrixCf[5];

    TH1F *hPtRecoNcCf[5];
    TH1F *hPtTrueMatchedNcCf[5];
    TH1F *hPtRecoMatchedNcCf[5];
    TH2F *hResponseMatrixNcCf[5];

    TH1F *hPtRecoTrackPlus[5];
    TH1F *hPtTrueMatchedTrackPlus[5];
    TH1F *hPtRecoMatchedTrackPlus[5];
    TH2F *hResponseMatrixTrackPlus[5];

    TH1F *hPtRecoTrackMinus[5];
    TH1F *hPtTrueMatchedTrackMinus[5];
    TH1F *hPtRecoMatchedTrackMinus[5];
    TH2F *hResponseMatrixTrackMinus[5];

    TH1F *hPtRecoClusterPlus[5];
    TH1F *hPtTrueMatchedClusterPlus[5];
    TH1F *hPtRecoMatchedClusterPlus[5];
    TH2F *hResponseMatrixClusterPlus[5];

    TH1F *hPtRecoClusterMinus[5];
    TH1F *hPtTrueMatchedClusterMinus[5];
    TH1F *hPtRecoMatchedClusterMinus[5];
    TH2F *hResponseMatrixClusterMinus[5];

    TH1F *hPtRecoTrClTight[5];
    TH1F *hPtTrueMatchedTrClTight[5];
    TH1F *hPtRecoMatchedTrClTight[5];
    TH2F *hResponseMatrixTrClTight[5];

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Trees
    //*******************************************************************************************************
    TTree *tTrueJets;
    std::vector<jets> true_jets;

    //*******************************************************************************************************

    TTree *tTotalParticles;
    std::vector<particles> total_particles;

    //*******************************************************************************************************

    TTree *tRecoJets;
    std::vector<jets> reco_jets;
    TTree *tMatchedJets;
    std::vector<jets> matched_jets;

    //*******************************************************************************************************

    inline void fillJESR(TH1F *hJESR[5][22], int cBin, float truePt, float factor)
    {
        if(truePt > PTBINS_TRUE_R3[0]  && truePt <= PTBINS_TRUE_R3[1])
            {
                hJESR[cBin][0]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[1]  && truePt <= PTBINS_TRUE_R3[2])
            {
                hJESR[cBin][1]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[2]  && truePt <= PTBINS_TRUE_R3[3])
            {
                hJESR[cBin][2]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[3]  && truePt <= PTBINS_TRUE_R3[4])
            {
                hJESR[cBin][3]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[4]  && truePt <= PTBINS_TRUE_R3[5])
            {
                hJESR[cBin][4]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[5]  && truePt <= PTBINS_TRUE_R3[6])
            {
                hJESR[cBin][5]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[6]  && truePt <= PTBINS_TRUE_R3[7])
            {
                hJESR[cBin][6]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[7]  && truePt <= PTBINS_TRUE_R3[8])
            {
                hJESR[cBin][7]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[8]  && truePt <= PTBINS_TRUE_R3[9])
            {
                hJESR[cBin][8]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[9]  && truePt <= PTBINS_TRUE_R3[10])
            {
                hJESR[cBin][9]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[10]  && truePt <= PTBINS_TRUE_R3[11])
            {
                hJESR[cBin][10]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[11]  && truePt <= PTBINS_TRUE_R3[12])
            {
                hJESR[cBin][11]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[12]  && truePt <= PTBINS_TRUE_R3[13])
            {
                hJESR[cBin][12]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[13]  && truePt <= PTBINS_TRUE_R3[14])
            {
                hJESR[cBin][13]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[14]  && truePt <= PTBINS_TRUE_R3[15])
            {
                hJESR[cBin][14]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[15]  && truePt <= PTBINS_TRUE_R3[16])
            {
                hJESR[cBin][15]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[16]  && truePt <= PTBINS_TRUE_R3[17])
            {
                hJESR[cBin][16]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[17]  && truePt <= PTBINS_TRUE_R3[18])
            {
                hJESR[cBin][17]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[18]  && truePt <= PTBINS_TRUE_R3[19])
            {
                hJESR[cBin][18]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[19]  && truePt <= PTBINS_TRUE_R3[20])
            {
                hJESR[cBin][19]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[20]  && truePt <= PTBINS_TRUE_R3[21])
            {
                hJESR[cBin][20]->Fill(factor);
            }
        if(truePt > PTBINS_TRUE_R3[21]  && truePt <= PTBINS_TRUE_R3[22])
            {
                hJESR[cBin][21]->Fill(factor);
            }
    }

};

#endif

























