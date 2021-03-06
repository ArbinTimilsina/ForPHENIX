#ifndef __CENTRALITYSTUDY_H__
#define __CENTRALITYSTUDY_H__

//Data classes I am using in analysis
#include <SubsysReco.h>
#include <PHCentralTrack.h>
#include <PHGlobal.h>
#include <EventHeader.h>
#include <PreviousEvent.h>

class PHCompositeNode;
class TH1F;
class TFile;

class CentralityStudy: public SubsysReco
{
 public:
  CentralityStudy();
  virtual ~CentralityStudy() {}

  //  For this analysis we only use following:
  int Init         (PHCompositeNode *topNode);
  int InitRun      (PHCompositeNode *topNode);
  int ResetEvent   (PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End          (PHCompositeNode *topNode);
  int EndRun       (const int runNumber);

 protected:
  TFile* outfile;
  unsigned int nTotalEvents;

  TH1F *hEvents;
  TH1F *hVertex;
  TH1F *hBbcQtotal;
};

#endif




