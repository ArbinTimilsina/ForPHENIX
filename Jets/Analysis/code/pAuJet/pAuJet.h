#ifndef __PAUJET_H__
#define __PAUJET_H__

//Data classes I am using in analysis
#include <SubsysReco.h>
#include <PHCentralTrack.h>
#include <PHGlobal.h>
#include <EventHeader.h>
#include <PreviousEvent.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>

//Cpp tools
#include <iomanip>

//Root tools
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCut.h>
#include <TVector3.h>
#include <TLorentzVector.h>

//For Trigger
#include <TrigLvl1.h>
#include <ErtOut.h>
#include <TriggerHelper.h>

//Used for EMCal warnmap
#include "EmcMap.h"

class PHCompositeNode;

#define VERTEX_CUT 30
#define TRACK_MIN_PT_CUT 0.5
#define CLUSTER_MIN_ENERGY_CUT 0.5

//pT bins
const int NPTBINS = 25;
const float PTBINS[NPTBINS + 1] =
    {
	5.0000, 5.4824, 6.0113, 6.5913, 7.2272, 7.9245, 8.6890, 9.5273, 10.4465, 11.4543, 12.5594, 13.7711,
	15.0998, 16.5566, 18.1539, 19.9054, 21.8258, 23.9315, 26.2404, 28.7720, 31.5479, 34.5915, 37.9289,
	41.5882, 45.6005, 50.0000
    };

namespace fastjet
{
    class JetDefinition;
}

typedef struct
{
    float mom;
    float theta;
    float pT;
    float eta;
    float phi;
    float phiDC;
    float zedDC;
    float alpha;
    float energy;
    float board;

    int charge;
    int quality;
    int n0;
    int arm;
    int armsect;
    int emcid;

    float px;
    float py;
    float pz;
    float eT;

    float pc3dphi;
    float pc3dz;
    float pc3sdphi;
    float pc3sdz;
    float emcdphi;
    float emcdz;
    float emcsdphi;
    float emcsdz;

} tracks;


typedef struct
{
    int arm;
    int sector;
    int armsect;
    int emcid;
    int yTowerPos;
    int zTowerPos;
    int towerId;

    bool ertFired;

    float energy;
    float theta;
    float pT;
    float eT;
    float phi;
    float px;
    float py;
    float pz;
    float eta;
    float chi2;
    float prob;

} clusters;

typedef struct
{
    int arm;

    float pT;
    float eta;
    float phi;
    float nc;
    float cf;

} jets;

class pAuJet: public SubsysReco
{
 public:
    pAuJet(std::string outfilename);
    virtual ~pAuJet() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void GetTracks   (PHCompositeNode *topNode, std::vector<tracks>& track_list);
    void GetClusters (PHCompositeNode *topNode, std::vector<clusters>& cluster_list, float vertex);

    void InitHistograms();

    void InitTrees(bool writeTrees);
    void FillTrees(bool writeTrees);

    float phiReduce(float phi)
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


 protected:
    const int verbo;
    const std::string outfname;
    const bool writeTree;

    TFile* outfile;

    PHGlobal *phglobal;
    EventHeader *evtheader;
    ErtOut *ertOut;
    TriggerHelper *triggerHelper;

    float zVertex;

    //For general information
    unsigned int runNumber;
    unsigned int eventNumber;

    unsigned int nTotalEvents;
    unsigned int nGoodEvents;

    unsigned int nJets;

    //List of tracks, clusters and jets
    TTree *tTracks;
    std::vector<tracks> tracks_list;

    TTree *tClusters;
    std::vector<clusters> clusters_list;

    TTree *tJets;
    std::vector<jets> jets_list;

    //Histograms
    TH1F *hVertex;
    TH1F *hEvents;

    TH1F *hPC3dphi;
    TH1F *hPC3dz;

    TH1F *hEMCdphi;
    TH1F *hEMCdz;

    TH2F *hTowerE;
    TH1F *hClustersTower;
    TH2F* hSectorHits[NARMSECT];

    TH1F *hCf;
    TH1F *hJets;
};
#endif

























