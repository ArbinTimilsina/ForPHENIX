#ifndef __EMCALMAPPP510_H__
#define __EMCALMAPPP510_H__

#include <SubsysReco.h>

class PHCompositeNode;
class TH3D;
class TFile;

#define NARMSECT 8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72

#define NECOREBIN 6

class EMCalMapPP510: public SubsysReco
{
 public:
  EMCalMapPP510(std::string _outfilename);
  virtual ~EMCalMapPP510() {}

  //  For this analysis we only use following:
  int Init         (PHCompositeNode *topNode);
  int InitRun      (PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End          (PHCompositeNode *topNode);
  int EndRun       (const int runNumber);
  int IsPbGl       (int armsect);

 protected:
  const int verbo;
  const std::string outfname;

  unsigned int runNumber;
  unsigned int nRunEvents;

  TH3D* hEmc3D[NARMSECT];

  TFile* outfile;
};

#endif
