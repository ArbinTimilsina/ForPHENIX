#ifndef __PC3EMCMATCHING_H__
#define __PC3EMCMATCHING_H__

#include "SubsysReco.h"
#include <TrackQuality.h>

//Used for PC3 matching
#include <Pc3Matching.h>

//Used for EMCal matching
#include <EmcMatching.h>

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
#define TRACK_MAX_PT_CUT 25.0
#define TRACK_MIN_PT_CUT 0.5

const int NARM = 2;

class Pc3EmcMatching: public SubsysReco
{
 public:
  Pc3EmcMatching(std::string _outfilename);
  virtual ~Pc3EmcMatching() {}

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

  TH2F *hEmcsdPhiFinal[NSECT][NZED][NCENT];
  TH2F *hEmcsdZFinal[NSECT][NZED][NCENT];

  TH2F *hPc3sdPhiFinal[NARM][NZED][NCENT];
  TH2F *hPc3sdZFinal[NARM][NZED][NCENT];

  TH1F *hAll;
  TH1F *hPc3;
  TH1F *hEmc;
  TH1F *hBoth;

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
