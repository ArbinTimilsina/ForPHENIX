#ifndef __ENERGYSCALECHECKS_H__
#define __ENERGYSCALECHECKS_H__

#include "JetAnalyzer.h"

//Data classes I am using in analysis
#include "Fun4AllServer.h"
#include "recoConsts.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1F;
class TFile;

class EnergyScaleChecks: public SubsysReco
{
 public:
    EnergyScaleChecks();
    virtual ~EnergyScaleChecks() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);
    void SetCuAu     (bool what);

    void GetPi0(std::vector<clusters> cluster_list,
                TH1F *hMass[8]);

 protected:
    bool isCuAu;

    JetAnalyzer *jetAnalyzer;
    TFile* outfile;

    unsigned int nTotalEvents;

    //*******************************************************************************************************
    //*******************************************************************************************************
    //Histograms
    //*******************************************************************************************************
    TH1F *hEvents;
    TH1F *hPi0Mass[8];
};

#endif




