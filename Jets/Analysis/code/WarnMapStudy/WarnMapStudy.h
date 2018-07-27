#ifndef __WARNMAPSTUDY_H__
#define __WARNMAPSTUDY_H__

#include "SubsysReco.h"
#include "Fun4AllHistoManager.h"
#include "EmcMap.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH2F;
class TH1F;
class TFile;

//Different cuts
#define VERTEX_CUT 10

class WarnMapStudy: public SubsysReco
{

 public:
  WarnMapStudy(std::string _outfilename);
  virtual ~WarnMapStudy() {}

  //  For this analysis we only use following:
  int Init         (PHCompositeNode *topNode);
  int InitRun      (PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End          (PHCompositeNode *topNode);
  int EndRun       (const int runNumber);

 protected:
  const int verbo;
  const std::string outfname;


  unsigned int runNumber;
  unsigned int nRunEvents;

  TH2F* hBefore[NARMSECT];
  TH2F* hAfter[NARMSECT];
  TH1F* hTower;

  TFile* outfile;
};

#endif /* __WARNMAPSTUDY_H__ */
