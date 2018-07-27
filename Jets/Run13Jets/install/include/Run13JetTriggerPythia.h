#ifndef __RUN13JETTRIGGERPYTHIA_H__
#define __RUN13JETTRIGGERPYTHIA_H__

//Data classes I am using in analysis
#include <SubsysReco.h>
#include <PHGlobal.h>
#include <EventHeader.h>
#include <PreviousEvent.h>

#include <TPythia6.h>
#include <TMath.h>
#include <TLorentzVector.h>

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
    int charge;

    float energy;
    float mom;
    float pT;
    float px;
    float py;
    float pz;
    float eta;
    float phi;

} pythiaParticles;

typedef struct
{
    float pT;
    float eta;
    float phi;

} pythiaJets;

class Run13JetTriggerPythia: public SubsysReco
{
 public:
    Run13JetTriggerPythia(const float iminPt, const int iTrigEvents);

    virtual ~Run13JetTriggerPythia() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);

    void GetPythia(PHCompositeNode *topNode,
                   std::vector<pythiaParticles>& particle_list,
                   std::vector<pythiaJets>& true_jet_list);

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
    float minPt;
    int nTrigEvents;

    unsigned int nTotalEvents;
    unsigned int nTriggeredEvents;
    unsigned int nTotalJets;
    unsigned int nTotalParticles;

    std::vector<pythiaParticles> pythia_particles;
    std::vector<pythiaJets> true_jets;
};

#endif
