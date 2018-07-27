#ifndef __JETNLOCORRECTION_H__
#define __JETNLOCORRECTION_H__

#include <JetAnalyzer.h>
#include <JetAnalyzerSim.h>

//Data classes I am using in analysis
#include <TPythia6.h>

//For Acceptance
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

    void GetPythia(PHCompositeNode *topNode,
                   std::vector<particles>& parton_list,
                   std::vector<particles>& final_list);

    void GetAntiKtParton(std::vector<particles> particle_list,
                         std::vector<jets>& antikt_jet_list);

    void GetAntiKtParticle(std::vector<particles> particle_list,
                           std::vector<jets>& antikt_jet_list);

    void GetMatchedJet(std::vector<jets> true_jet_list,
                       std::vector<jets> good_jet_list,
                       std::vector<jets>& matched_jet_list,
                       TH1F *hDr,
                       TH1F *hMatched,
                       TH1F *hEfficiency,
                       TH2F *hMatrix,
                       TH1F *hForCF[17]);
 protected:
    const std::string outfname;

    TFile* outfile;

    JetAnalyzer *jetAnalyzer;

    unsigned int nRunEvents;

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    //General
    //*****************************************************************************************************************************
    TH1F *hEvents;

    TH1F *hPtTrueJet;
    TH1F *hPtRecoJet;

    TH1F *hDistance;
    TH1F *hPtMatchedJet;
    TH1F *hPtForEfficiency;
    TH2F *hResponseMatrix;

    TH1F *hForCorrectionFactor[17];

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Ntuples
    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Pythia
    TTree *pythiaParton;
    std::vector<particles> pythia_parton_list;

    TTree *pythiaFinalParticle;
    std::vector<particles> pythia_final_particle_list;

    //Jets
    TTree *pythiaPartonJet;
    std::vector<jets> pythia_parton_jet_list;

    TTree *pythiaFinalParticleJet;
    std::vector<jets> pythia_final_particle_jet_list;

    TTree *pythiaMatchedJet;
    std::vector<jets> pythia_matched_jet_list;

};

#endif













