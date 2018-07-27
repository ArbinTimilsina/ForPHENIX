#ifndef __FOREMCALMATCHING_H__
#define __FOREMCALMATCHING_H__

#include "SubsysReco.h"
#include "TrackQuality.h"

//Used for PC3 mathing
#include "CuAu200_PC3_matching.h"

//Data classes I am using in analysis
#include "PHCentralTrack.h"
#include "PHGlobal.h"
#include "EventHeader.h"
#include "PreviousEvent.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "TriggerHelper.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1F;
class TH2F;
class TFile;
class TNtuple;

#define VERTEX_CUT 10

class forEMCalMatching: public SubsysReco
{
 public:
  forEMCalMatching(std::string _outfilename);
  virtual ~forEMCalMatching() {}

  //  For this analysis we only use following:
  int Init         (PHCompositeNode *topNode);
  int InitRun      (PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent   (PHCompositeNode *topNode);
  int End          (PHCompositeNode *topNode);
  int EndRun       (const int runNumber);

 protected:
  const int verbo;
  const std::string outfname;

  TFile* outfile;

  TNtuple* goodTracks;

  //For centrality bin
  unsigned int cbin;

  //For general information
  unsigned int runNumber;
  unsigned int nEvents;
  unsigned int eventNumber;

  PHGlobal *phglobal;
  EventHeader *evtheader;
  PHCentralTrack *phcentral;

  float zvertex;
  float bbct0;
  float centrality;

};

#endif
