#ifndef __EMCALMAP_H__
#define __EMCALMAP_H__

#include <SubsysReco.h>

class PHCompositeNode;
class TH3D;
class TFile;

//Vertex cut
#define VERTEX_CUT 10

#define NARMSECT 8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72

#define NECOREBIN 8

class EMCalMap: public SubsysReco
{
 public:
  EMCalMap(std::string _outfilename);
  virtual ~EMCalMap() {}

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

#endif /* __EMCALMAP_H__ */
