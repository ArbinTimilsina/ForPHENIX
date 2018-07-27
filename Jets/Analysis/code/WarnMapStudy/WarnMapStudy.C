//General PHENIX tools
#include "getClass.h"
#include "PHCompositeNode.h"
#include "phool.h"
#include "RunHeader.h"

//Fun4All tools
#include "Fun4AllServer.h"
#include "Fun4AllHistoManager.h"
#include "Fun4AllReturnCodes.h"

//C tools
#include <cstdlib>

//  Root tools
#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include "TVector3.h"
#include "TLorentzVector.h"

//  Data classes I am using in analysis
#include "PHCentralTrack.h"
#include "PHGlobal.h"
#include "PreviousEvent.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"

//My source file
#include "WarnMapStudy.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables

WarnMapStudy::WarnMapStudy(string _outfilename)
  : SubsysReco("WarnMapStudy"),
    verbo(1),
    outfname(_outfilename)
{
  return;
}

int WarnMapStudy::Init(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  outfile = new TFile(outfname.c_str(), "RECREATE");

  //load and read the warnmap
  ReadWarnMap();

  for (int ias = 0; ias < NARMSECT; ias++)
    {
      char hnbefore[256];
      sprintf(hnbefore, "hBefore_sector_%i", ias);

      char htbefore[256];
      sprintf(htbefore, "Before secotr %i", ias);

      char hnafter[256];
      sprintf(hnafter, "hAfter_sector_%i", ias);

      char htafter[256];
      sprintf(htafter, "After sector %i", ias);

      int ny, nz;
      if (IsPbGl(ias))
        {
	  ny = YPOS_PBGL;
          nz = ZPOS_PBGL;
        }
      else
        {
          ny = YPOS_PBSC;
          nz = ZPOS_PBSC;
        }

      hBefore[ias] = new TH2F(hnbefore, htbefore, nz, 0, nz, ny, 0, ny);
      hAfter[ias] = new TH2F(hnafter, htafter, nz, 0, nz, ny, 0, ny);
    }

  hTower = new TH1F("hTower","Hit vs Tower-id",NTOWER,0,NTOWER);

  return 0;
}

int WarnMapStudy::InitRun(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  nRunEvents = 0;

  RunHeader *run_header = getClass<RunHeader> (topNode, "RunHeader");
  if (!run_header)
    {
      cout << "No RunHeader! No sense continuing" << endl;
      exit(1);
    }
  runNumber = run_header->get_RunNumber();

  return EVENT_OK;
}

int WarnMapStudy::process_event(PHCompositeNode *topNode)
{
  // Informational message...
  nRunEvents++;
  if (nRunEvents % 1000 == 0 && verbosity)
    {
      if (verbo) cout << "Events for run " << runNumber << " = " << nRunEvents << endl;
    }

  //  Get the data I need...
  PHGlobal *phglobal                     = getClass<PHGlobal>                (topNode, "PHGlobal");
  if (nRunEvents == 1 && !phglobal)
    {
      cout << "No PHGlobal!  No sense continuing" << endl;
      exit(1);
    }

  PHCentralTrack *phcentral              = getClass<PHCentralTrack>          (topNode, "PHCentralTrack");
  if (nRunEvents == 1 && !phcentral)
    {
      cout << "No PHCentral!  No sense continuing" << endl;
      exit(1);
    }

  emcClusterContainer *emcclustercontainer = getClass<emcClusterContainer>      (topNode, "emcClusterContainer");
  if (nRunEvents == 1 && !emcclustercontainer)
    {
      cout << "No emcClusterContainer!  No sense continuing" << endl;
      exit(1);
    }

  double zvertex = phglobal->getBbcZVertex();

  //Vertex cut
  if (fabs(zvertex) > VERTEX_CUT) return DISCARDEVENT;


  //The number of towers is 48(y) x 96(z) x 2(PbGl) + 36(y) x 72(z) x 6(PbSc) = 9216+15552 = 24768.
  //Tower ID is a unique ID of all towers in eight sectors (0~24767).

  int Nclus = emcclustercontainer->size();
  for (int iclus = 0; iclus < Nclus; iclus++)
    {
      emcClusterContent* clus = emcclustercontainer->getCluster(iclus);

      int arm  = clus->arm(); //In EMCal convention, West Arm is 0 and East Arm is 1, and thus armsector 0...7 are W0...W3 E0...E3
      int sector = clus->sector();
      int armsects = (arm * 4) + sector;
      int yPos = clus->iypos();
      int zPos = clus->izpos();
      float ecore = clus->ecore();

      //Ecore cut
      if(ecore<0.6) continue;
      hBefore[armsects]->Fill(zPos, yPos,1);

      if(IsHot(armsects,yPos,zPos) || IsDead(armsects,yPos,zPos) || IsUncalib(armsects,yPos,zPos)) continue;
      hAfter[armsects]->Fill(zPos, yPos,1);

      int towerId = GetTowerID(armsects,yPos,zPos);
      hTower->Fill(towerId,1);
    }

  // any other return code might lead to aborting the event or analysis
  return 0;
}

int WarnMapStudy::EndRun(const int runNumber)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  if (verbo)
    {
      cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
      cout << "Run number:                                 " << runNumber << endl;
      cout << "Number of total events for this run:        " << nRunEvents << endl;

    }
  return 0;
}

int WarnMapStudy::End(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  End called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  outfile->Write(outfname.c_str());
  outfile->Close();

  if (verbo)
    {
      cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
    }
  return 0;
}

