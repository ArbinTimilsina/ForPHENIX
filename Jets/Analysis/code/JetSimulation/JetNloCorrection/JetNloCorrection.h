#ifndef __JETNLOCORRECTION_H__
#define __JETNLOCORRECTION_H__

#include <JetTriggerPythia.h>

//Data classes I am using in analysis
#include <TPythia6.h>
#include "TAcceptParticle.h"

class PHCompositeNode;
class PHPythiaContainer;
class TH1F;
class TFile;
class TTree;

namespace fastjet
{
    class JetDefinition;
}

class JetNloCorrection: public SubsysReco
{
 public:
    JetNloCorrection();
    virtual ~JetNloCorrection() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);

    void GetParton(std::vector<particles>& parton_list);

    void GetAntiKtParton(std::vector<particles> particle_list,
                         std::vector<jets>& antikt_jet_list);

    void GetAntiKtParticle(PHCompositeNode *topNode,
                           std::vector<particles>& particle_list,
                           std::vector<jets>& antikt_jet_list);

 protected:
    const std::string outfname;

    TFile* outfile;

    PHPythiaContainer *phpythia;

    unsigned int nRunEvents;
    bool writeTree;

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    //General
    //*****************************************************************************************************************************
    TH1F *hEvents;

    TH1F *hPtParton;
    TH1F *hPtParticle;

    TTree *pythiaParton;
    std::vector<particles> pythia_parton_list;

    TTree *pythiaParticle;
    std::vector<particles> pythia_particle_list;

    TTree *pythiaPartonJet;
    std::vector<jets> pythia_parton_jet_list;

    TTree *pythiaParticleJet;
    std::vector<jets> pythia_particle_jet_list;

    inline float phiReduce(float phi)
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

    inline bool passPhenixAcceptance(float jetEta, float jetPhi)
    {
        return((fabs(jetEta) < 0.35) && ((jetPhi > -0.589 && jetPhi < 0.982) || (jetPhi > 2.160 && jetPhi < 3.731)));
    }
};

#endif
















