#ifndef __JETHADCORRECTION_H__
#define __JETHADCORRECTION_H__

//Data classes I am using in analysis
#include <SubsysReco.h>
#include <TPythia6.h>
#include <TMath.h>

class PHCompositeNode;
class PHPythiaContainer;
class TH1F;
class TFile;
class TTree;

const int NPTBINS = 16;
const float PTBINS[NPTBINS + 1] =
    {
	9.5273,	10.4465, 11.4543, 12.5594, 13.7711, 15.0998, 16.5566, 18.1539, 19.9054,
	21.8258, 23.9315, 26.2404, 28.7720, 31.5479, 34.5915, 37.9289, 41.5882
    };

namespace fastjet
{
    class JetDefinition;
}

typedef struct
{
    int arm;
    int charge;
    float energy;
    float mom;
    float pT;
    float eT;
    float px;
    float py;
    float pz;
    float eta;
    float phi;
} particles;

class JetHadCorrection: public SubsysReco
{
 public:
    JetHadCorrection();
    virtual ~JetHadCorrection() {}

    //  For this analysis we only use following:
    int Init         (PHCompositeNode *topNode);
    int InitRun      (PHCompositeNode *topNode);
    int ResetEvent   (PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End          (PHCompositeNode *topNode);

    void GetAntiKtJets(std::vector<particles> particle_list, float R, TH1F *hJet);


 protected:
    const std::string outfname;

    TFile* outfile;
    PHPythiaContainer *phpythia;

    unsigned int nRunEvents;
    bool writeTree;

    TH1F *hEvents;

    TH1F *hPtPartonJetsR2;
    TH1F *hPtParticleJetsR2;

    TH1F *hPtPartonJetsR3;
    TH1F *hPtParticleJetsR3;

    TTree *pythiaParton;
    std::vector<particles> pythia_parton_list;

    TTree *pythiaParticle;
    std::vector<particles> pythia_particle_list;

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
};

#endif
















