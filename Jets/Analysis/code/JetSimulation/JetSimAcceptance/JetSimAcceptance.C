//General PHENIX tools
#include "PHCompositeNode.h"
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include "phool.h"
#include <getClass.h>
#include "RunHeader.h"

//Fun4All tools
#include "Fun4AllServer.h"
#include "Fun4AllHistoManager.h"
#include "Fun4AllReturnCodes.h"

//Root tools
#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include <TNtuple.h>

#include <VtxOut.h>

//My source file
#include "JetSimAcceptance.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
JetSimAcceptance::JetSimAcceptance()
  : SubsysReco("JetSimAcceptance")
{
  nTotalEvents = 0;

  return;
}

int JetSimAcceptance::Init(PHCompositeNode *topNode)
{
  //For data code
  jetAnalyzer = new JetAnalyzer("dummyFile.root");
  jetAnalyzer->SetData(true);
  if(isCuAu)
    {
      jetAnalyzer->SetCuAu(true);
    }
  else
    {
      jetAnalyzer->SetCuAu(false);
    }
  jetAnalyzer->MyInit();

  //Output file name
  outfile = new TFile("JetSimAcceptance.root", "RECREATE");

  //*******************************************************************************************************
  //*******************************************************************************************************
  //Histograms
  //*******************************************************************************************************
  hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

  //*******************************************************************************************************
  hModifiedQuality_NE = new TH2F("hModifiedQuality_NE", "Modified Quality, NE",  400, 0, 80, 120, -0.6, 0.6);
  hModifiedQuality_SE = new TH2F("hModifiedQuality_SE", "Modified Quality, SE",  400, 0, 80, 120, -0.6, 0.6);
  hModifiedQuality_NW = new TH2F("hModifiedQuality_NW", "Modified Quality, NW",  400, 0, 80, 120, -0.6, 0.6);
  hModifiedQuality_SW = new TH2F("hModifiedQuality_SW", "Modified Quality, SW",  400, 0, 80, 120, -0.6, 0.6);

  hTracksPhi = new TH1F("hTracksPhi", "#frac{dN}{d#phi} of tracks going into Jet Reconstruction", 200, -1.0, 4.0);

  //*******************************************************************************************************
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      int ny, nz;
      if (EmcMap::IsPbGl(ias))
        {
	  ny = YPOS_PBGL;
	  nz = ZPOS_PBGL;
        }
      else
        {
	  ny = YPOS_PBSC;
	  nz = ZPOS_PBSC;
        }
      hSectorHits[ias] = new TH2F(Form("hSecorHits_%i", ias), Form("Sector %i", ias), nz, 0, nz, ny, 0, ny);
    }

  hClustersPhi = new TH1F("hClustersPhi", "#frac{dN}{d#phi} of clusters going into Jet Reconstruction", 200, -1.0, 4.0);

  //*******************************************************************************************************
  return EVENT_OK;
}

int JetSimAcceptance::InitRun(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int JetSimAcceptance::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}


int JetSimAcceptance::process_event(PHCompositeNode *topNode)
{
  nTotalEvents++;

  getAcceptance(topNode, hModifiedQuality_NE, hModifiedQuality_SE, hModifiedQuality_NW, hModifiedQuality_SW, hTracksPhi,
		hSectorHits, hClustersPhi);
  return 0;
}

void JetSimAcceptance::SetCuAu(bool what = true)
{
  isCuAu = what;
}

void JetSimAcceptance::getAcceptance(PHCompositeNode *topNode,
                                     TH2F *hNE, TH2F *hSE, TH2F *hNW, TH2F *hSW, TH1F *hTrackPhi,
                                     TH2F *hHits[NARMSECT], TH1F *hClusterPhi)

{
    VtxOut *vertexOut = getClass<VtxOut>(topNode, "VtxOut");
    if (!vertexOut)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    float zvertex = vertexOut->get_ZVertex();

  //Tracks
  std::vector<tracks> track_list;
  track_list.clear();
  jetAnalyzer->GetTracks(topNode, track_list);

  for (unsigned int t = 0; t < track_list.size(); t++)
    {
      float board        = track_list[t].board;
      bool passQuality   = track_list[t].passQuality;
      float pT           = track_list[t].pT;
      float alpha        = track_list[t].alpha;
      float zed          = track_list[t].zed;
      int arm            = track_list[t].arm;

      if (pT > 0.2 && pT < 25)
        {
	  if (passQuality)
            {
	      if (zed > 0 && arm == 0)
                {
		  hNE->Fill(board, alpha);
                }
	      if (zed > 0 && arm == 1)
                {
		  hNW->Fill(board, alpha);
                }
	      if (zed < 0 && arm == 0)
                {
		  hSE->Fill(board, alpha);
                }
	      if (zed < 0 && arm == 1)
                {
		  hSW->Fill(board, alpha);
                }
            }
        }
    }

  //Cluster
  std::vector<clusters> cluster_list;
  cluster_list.clear();
  jetAnalyzer->SetData(false);
  jetAnalyzer->GetClusters(topNode, cluster_list, zvertex);

  for (int iclus = 0; iclus < cluster_list.size(); iclus++)
    {
      if(cluster_list[iclus].passEverything)
        {
	  int armsect = cluster_list[iclus].armsect;
	  int yPos = cluster_list[iclus].yPos;
	  int zPos = cluster_list[iclus].zPos;
	  hHits[armsect]->Fill(zPos, yPos);
        }
    }


  std::vector<tracks> charged_particle_list;
  charged_particle_list.clear();

  std::vector<clusters> neutral_particle_list;
  neutral_particle_list.clear();

  std::vector<particles> particle_list;
  particle_list.clear();

  jetAnalyzer->GetParticles(track_list, cluster_list, charged_particle_list, neutral_particle_list, particle_list);

  for(unsigned int c = 0; c < charged_particle_list.size(); c++)
    {
      hTrackPhi->Fill(charged_particle_list[c].phi);
    }

  for(unsigned int n = 0; n < neutral_particle_list.size(); n++)
    {
        hClusterPhi->Fill(neutral_particle_list[n].phi);
    }
}


int JetSimAcceptance::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetSimAcceptance::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nTotalEvents);

    cout << endl;
    cout << "Total events processed: " << nTotalEvents << endl;
    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}

