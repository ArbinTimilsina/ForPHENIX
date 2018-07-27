#ifndef __RUNQA_H__
#define __RUNQA_H__

#include "SubsysReco.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH2F;
class TH1F;
class TFile;
class TNtuple;
class TriggerHelper;

#define IRUN_NUMBER 357500
#define FRUN_NUMBER 377400
#define RUN_BIN 19900
#define NARMSECTS 8

class RunQA: public SubsysReco
{

 public:
    RunQA(std::string _outfilename);
    virtual ~RunQA() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

 protected:
    const int verbo;
    const std::string outfname;

    TFile* outfile;

    unsigned int runNumber;
    unsigned int nRunEvents;

    double runVertex;
    double runCentrality;

    unsigned int nTracksT;
    unsigned int nTracksQ;
    unsigned int nTracksE;
    unsigned int nTracksW;
    unsigned int nEMCalHits[NARMSECTS];

    //Run QA histograms
    TH1F *hEvents;

    TH1F *hVertex;
    TH1F *hCentrality;

    TH1F *hTracksT;
    TH1F *hTracksQ;
    TH1F *hTracksE;
    TH1F *hTracksW;

    TH1F *hEMCalHits[NARMSECTS];
};

#endif





