#ifndef __JETSIMACCEPTANCE_H__
#define __JETSIMACCEPTANCE_H__

#include "JetAnalyzer.h"

//Data classes I am using in analysis
#include "Fun4AllServer.h"
#include "recoConsts.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1F;
class TFile;

class JetSimAcceptance: public SubsysReco
{
 public:
  JetSimAcceptance();
  virtual ~JetSimAcceptance() {}

  //  For this analysis we only use following:
  int Init         (PHCompositeNode *topNode);
  int InitRun      (PHCompositeNode *topNode);
  int ResetEvent   (PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End          (PHCompositeNode *topNode);
  int EndRun       (const int runNumber);
  void SetCuAu     (bool what);

  void getAcceptance(PHCompositeNode *topNode,
		     TH2F *hNE, TH2F *hSE, TH2F *hNW, TH2F *hSW, TH1F *hTrackPhi,
		     TH2F *hHits[NARMSECT], TH1F *hClusterPhi);

 protected:
  bool isCuAu;

  TFile* outfile;

  unsigned int nTotalEvents;

  JetAnalyzer *jetAnalyzer;

  //*******************************************************************************************************
  //*******************************************************************************************************
  //Histograms
  //*******************************************************************************************************
  TH1F *hEvents;

  //*******************************************************************************************************
  TH2F *hModifiedQuality_NE;
  TH2F *hModifiedQuality_SE;
  TH2F *hModifiedQuality_NW;
  TH2F *hModifiedQuality_SW;

  TH1F *hTracksPhi;
  //*******************************************************************************************************

  TH2F* hSectorHits[NARMSECT];

    TH1F *hClustersPhi;
};

#endif




