#ifndef __ALPHAVSPHI_H__
#define __ALPHAVSPHI_H__

#include <string>
#include <SubsysReco.h>
#include <TFile.h>

class PHCompositeNode;
class PHCentralTrack;
class PHGlobal;
class TH2F;

//For momentum
const float MomfactorWest = 0.9742;
const float MomfactorEast = 0.9727;

class AlphaVsPhi : public SubsysReco
{
 public:
    AlphaVsPhi(std::string outfilename);
    virtual ~AlphaVsPhi() {}

    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SetCuAu(bool what);

 protected:
    const std::string outfname;
    TFile* outfile;

    unsigned int nEvents;
    bool isCuAu;

    TH2F *hAlphaVsPhi_East;
    TH2F *hAlphaVsPhi_West;

    inline float getNewAlpha(float alpha, float phi, float dcarm)
    {
        float xp = sin(phi);
        float yp = cos(phi);

        if (dcarm == 0)
            {
                float XOffsetE = 1.19326e-01;
                float YOffsetE = 1.30366e-01;

                float AlphaOffsetE = (XOffsetE * xp / 220. + YOffsetE * yp / 220.);
                return (alpha - AlphaOffsetE);
            }
        else if (dcarm == 1)
            {
                float XOffsetW = 1.81180e-01;
                float YOffsetW = 1.48351e-01;

                float AlphaOffsetW = (XOffsetW * xp / 220. + YOffsetW * yp / 220.);
                return (alpha - AlphaOffsetW);
            }
        return alpha;
    }

};
#endif
