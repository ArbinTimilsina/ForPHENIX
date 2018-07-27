#ifndef __JETCONTAINERANALYZER_H__
#define __JETCONTAINERANALYZER_H__

#include <JetTriggerPythia.h>
#include <JetAnalyzer.h>

#include <PHHijing.h>
#include <PHHijingHeader.h>
#include <PHHijingHeaderV3.h>
#include <PHPythia.h>
#include <PHPythiaContainer.h>
#include <PHPythiaContainerV2.h>

//Data classes I am using in analysis
#include <Fun4AllServer.h>
#include <recoConsts.h>

const float MM2CM = 0.1;

//class Fun4AllHistoManager;
class PHCompositeNode;
class TH1F;
class TFile;
class TTree;

namespace fastjet
{
    class JetDefinition;
}

class JetContainerAnalyzer: public SubsysReco
{
 public:
    JetContainerAnalyzer(const float ir, const float iminPt,
                         const float incPythia, const float iminCfPythia, const float imaxCfPythia);

    virtual ~JetContainerAnalyzer() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);
    int EndRun       (const int runNumber);
    void SetCuAu     (bool what);

    void getPythiaTrueJet(std::vector<particles>& particle_list,
                          std::vector<jets>& true_jet_list);

    void getPisaRecoParticles(std::vector<particles>& good_particle_list);

    void getTotalParticles(std::vector<particles> sim_good_particle_list,
                           std::vector<particles>& total_particle_list);


 protected:
    int CreateNodeTree(PHCompositeNode *topNode);

    float R;
    float minPt;
    float ncPythia;
    float minCfPythia;
    float maxCfPythia;

    bool isCuAu;

    Fun4AllServer* se;
    recoConsts* rc;

    JetAnalyzer *jetAnalyzer;
    PHPythiaContainer *phpythia;
    JetTriggerPythia *jetTriggerPythia;
    PHGlobal *phglobal;

    PHCompositeNode *dataNode;
    PHCompositeNode *pythiaNode;
    PHCompositeNode *pisaRecoNode;

    unsigned int nTotalEvents;
    unsigned int nGoodEvents;

    float centrality;
    float zVertexPythia;
    float zVertexPisaReco;
    float zVertexData;

    PHHijingHeader *eventContainer;
    PHPythiaContainer *jetContainer;
    PHPythiaContainer *particleContainer;

    //*******************************************************************************************************
    std::vector<particles> pythia_particles;
    std::vector<jets> true_jets;
    std::vector<particles> pisa_reco_good_particles;
    std::vector<particles> total_particles;
};

#endif

























