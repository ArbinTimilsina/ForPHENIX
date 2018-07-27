#ifndef __JETTRIGGERPYTHIA_H__
#define __JETTRIGGERPYTHIA_H__

#include <JetAnalyzer.h>

//Data classes I am using in analysis
#include <TPythia6.h>

//For Acceptance
#include "TAcceptParticle.h"

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


class JetTriggerPythia: public SubsysReco
{
 public:
    JetTriggerPythia(const float ir, const float inc, const float iminPt, const float iminCf, const float imaxCf,
                     const int iTrigEvents);

    virtual ~JetTriggerPythia() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void GetPythia(PHCompositeNode *topNode,
                   std::vector<particles>& particle_list,
                   std::vector<jets>& true_jet_list);
 protected:
    float R;
    float nc;
    float minPt;
    float minCf;
    float maxCf;

    int nTrigEvents;

    unsigned int nTotalEvents;
    unsigned int nTriggeredEvents;
    unsigned int nTotalJets;

    unsigned int nTotalParticles;

    JetAnalyzer *jetAnalyzer;

    std::vector<particles> pythia_particles;
    std::vector<jets> true_jets;
};

#endif











