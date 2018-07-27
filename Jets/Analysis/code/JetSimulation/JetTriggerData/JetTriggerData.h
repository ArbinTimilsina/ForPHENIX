#ifndef __JETTRIGGERDATA_H__
#define __JETTRIGGERDATA_H__

#include <JetAnalyzer.h>
#include <iostream>
#include <fstream>

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

class JetTriggerData: public SubsysReco
{
 public:
    JetTriggerData(const float ir, const float inc, const float iminPt, const float iminCf, const float imaxCf,
                   const int iTrigEvents);
    virtual ~JetTriggerData() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

 protected:
    float R;
    float nc;
    float minPt;
    float minCf;
    float maxCf;
    int nTrigEvents;

    float zvertex;
    float centrality;

    ofstream vertexFile;

    unsigned int nTotalEvents;
    unsigned int nTotalFailedVertexEvents;
    unsigned int nTotalJetEvents;
    unsigned int nTriggeredEvents;

    JetAnalyzer *jetAnalyzer;
};

#endif











