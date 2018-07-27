//General PHENIX tools
#include "getClass.h"
#include "PHCompositeNode.h"
#include "phool.h"
#include "RunHeader.h"

//Fun4All tools
#include "Fun4AllServer.h"
#include "Fun4AllHistoManager.h"
#include "Fun4AllReturnCodes.h"

//  Root tools
#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TNtuple.h>

//My source file
#include "forEMCalMatching.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables

forEMCalMatching::forEMCalMatching(string _outfilename)
  : SubsysReco("forEMCalMatching"),
    verbo(1),
    outfname(_outfilename)
{
  return;
}

int forEMCalMatching::Init(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  outfile = new TFile(outfname.c_str(), "RECREATE");

  //For PC3 matching
  pc3_init_fit_pars_II();

  string tracks = "armsect:zed:centrality:pT:charge:phi:emcdphi:emcdz:theta:beta";
  goodTracks = new TNtuple("goodTracks", "Good Tracks for EMC matching", tracks.c_str());

  return 0;
}

int forEMCalMatching::InitRun(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>> InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  nEvents = 0;

  RunHeader *run_header = getClass<RunHeader> (topNode, "RunHeader");
  if (!run_header)
    {
      cout << "No RunHeader! No sense continuing" << endl;
      exit(1);
    }

  runNumber = run_header->get_RunNumber();

  return EVENT_OK;
}

int forEMCalMatching::process_event(PHCompositeNode *topNode)
{
  // Informational message...
  nEvents++;
  if (nEvents % 1000 == 0 && verbosity)
    {
      if (verbo) cout << "Events for run " << runNumber << " = " << nEvents << endl;
    }

  //  Get the data I need...
  phglobal                               = getClass<PHGlobal>                (topNode, "PHGlobal");
  if (nEvents == 1 && !phglobal)
    {
      cout << "No PHGlobal!  No sense continuing" << endl;
      exit(1);
    }

  evtheader                              = getClass<EventHeader>             (topNode, "EventHeader");
  if (nEvents == 1 && !evtheader)
    {
      cout << "No EventHeader! No sense continuing" << endl;
      exit(1);
    }

  phcentral                              = getClass<PHCentralTrack>          (topNode, "PHCentralTrack");
  if (nEvents == 1 && !phcentral)
    {
      cout << "No PHCentral!  No sense continuing" << endl;
      exit(1);
    }

  zvertex = phglobal->getBbcZVertex();

  bbct0 = phglobal->getBbcTimeZero();

  centrality = phglobal->getCentrality();

  //Vertex cut
  if (fabs(zvertex) > VERTEX_CUT) return DISCARDEVENT;

  //Centrality cut
  if (centrality < 0 || centrality > 100) return DISCARDEVENT;

  if (evtheader) eventNumber = evtheader->get_EvtSequence();

  for (unsigned int i = 0; i < phcentral->get_npart(); i++)
    {
      float mom           = phcentral->get_mom(i); // Magnitude of the momentum.
      float theta0        = phcentral->get_the0(i); //The track's theta direction at the vertex
      float pT            = mom * sin(theta0);
      int quality         = phcentral->get_quality(i); // Quality of the Drift Chamber Tracks
      int arm             = phcentral->get_dcarm(i); //Arm containing the track (East=0, West=1)
      int n0              = phcentral->get_n0(i); //The number of phototubes that fired in the normally sized ring area
      float zed           = phcentral->get_zed(i); //The Z coordinate at which the track crosses PC1
      int emcid           = phcentral->get_emcid(i); // Index of the emc cluster used
      float ecore         = phcentral->get_ecore(i); // EMC "shower core" energy
      float phi           = phcentral->get_phi(i); // The phi coordinate at which the track crosses the drift chamber reference radius
      float phi0          = phcentral->get_phi0(i); // The track's phi direction at the vertex
      float alpha         = phcentral->get_alpha(i); // This is the inclination of the track w.r.t. an infinite momentum track
      float pc3dphi       = phcentral->get_pc3dphi(i);
      float pc3dz         = phcentral->get_pc3dz(i);
      float emcdphi       = phcentral->get_emcdphi(i);
      float emcdz         = phcentral->get_emcdz(i);
      float eta           = (float)(-log(tan(0.5 * theta0)));
      float charge        = phcentral->get_charge(i);
      int sector          = phcentral->get_sect(i); //EMC sector of the associated cluster
      float beta          = phcentral->get_beta(i);

      // Charge not found
      bool passCharged = !((charge == -9999) || (charge == 0) || (alpha == -9999) || (alpha == 0));

      //DC zed cut
      bool passZed =  (fabs(zed) < 75.0 );

      bool passMost = passCharged && passZed;

      if (!passMost) continue;

      //For armsect:
      //  3 (W3PbSc) || 7 (E3PbSc)
      //  2 (W2PbSc) || 6 (E2PbSc)
      //  1 (W1PbSc) || 5 (E1PbGl)
      //  0 (W0PbSc) || 4 (E0PbGl)

      int armsect = -99;
      if (arm == 1)
        {
          armsect = sector;
        }
      else
        {
          armsect = 4 + sector;
        }
      if (armsect < 0 || armsect >= 8) continue;

      if (emcdphi == -9999 || emcdz == -9999) continue;

      float board = -9999;
      if (arm == 0) board = (3.72402 - phi + 0.008047 * cos(phi + 0.87851)) / 0.01963496;
      if (arm == 1) board = (0.573231 + phi + 0.0046 * cos(phi + 0.05721)) / 0.01963496;

      //Quality cut
      bool passQuality1 = (quality == 63) || (quality == 31);
      bool passQuality2 = TrackQuality::passQualityMask(quality, board, alpha, zed, arm);

      //PC3 Matching
      float pc3sdphi = -99.9;
      float pc3sdz   = -99.9;

      int ized = -9999;
      if     (zed > -70 && zed <= -57) ized = 0;
      else if (zed > -57 && zed <= -44) ized = 1;
      else if (zed > -44 && zed <= -31) ized = 2;
      else if (zed > -31 && zed <= -18) ized = 3;
      else if (zed > -18 && zed <= -5) ized = 4;
      else if (zed >= 5  && zed <  18) ized = 5;
      else if (zed >= 18 && zed <  31) ized = 6;
      else if (zed >= 31 && zed <  44) ized = 7;
      else if (zed >= 44 && zed <  57) ized = 8;
      else if (zed >= 57 && zed <  70) ized = 9;

      if (!(ized < 0 || ized > 9))
        {
          int NMUL = 10;
          int icent = (int) ( NMUL * ((centrality - 0.001) / 100.) );
          if (!(icent < 0 || icent > NMUL))
            {
              pc3sdphi = pc3_sdphi_func_II(charge, arm, icent, ized, pT, pc3dphi);
              pc3sdz   = pc3_sdz_func_II(charge, arm, icent, ized, pT, pc3dz);
            }
        }

      bool pc3Matching = sqrt((pc3sdphi * pc3sdphi) + (pc3sdz * pc3sdz)) < 3.0;

      float ntdata[10];

      ntdata[0] = armsect;
      ntdata[1] = zed;
      ntdata[2] = centrality;
      ntdata[3] = pT;
      ntdata[4] = charge;
      ntdata[5] = phi;
      ntdata[6] = emcdphi;
      ntdata[7] = emcdz;
      ntdata[8] = theta0;
      ntdata[9] = beta;

      if ((pT > 0.2 && pT <= 3.0 && passQuality1) || (pT > 3.0 && pT<25.0 && passQuality2 && pc3Matching && n0>0))
        {
          goodTracks->Fill(ntdata);
        }
    }

  // any other return code might lead to aborting the event or analysis
  return 0;
}

int forEMCalMatching::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int forEMCalMatching::EndRun(const int runNumber)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
  return EVENT_OK;
}

int forEMCalMatching::End(PHCompositeNode *topNode)
{
  if (verbo) cout << ">>>>>>>>>>>>>>>>>>>  End called <<<<<<<<<<<<<<<<<<<<<<<" << endl;

  outfile->Write(outfname.c_str());
  outfile->Close();

  if (verbo)
    {
      cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
    }
  return EVENT_OK;
}

