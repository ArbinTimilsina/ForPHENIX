//General PHENIX tools
#include <getClass.h>
#include <PHCompositeNode.h>
#include <phool.h>
#include <RunHeader.h>
#include <PHGlobal.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//C tools
#include <cstdlib>

//  Root tools
#include <TH3D.h>
#include <TFile.h>

//My source file
#include "EMCalMapPP510.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables

EMCalMapPP510::EMCalMapPP510(string _outfilename)
    : SubsysReco("EMCalMapPP510"),
      verbo(1),
      outfname(_outfilename)
{
    return;
}

int EMCalMapPP510::Init(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    outfile = new TFile(outfname.c_str(), "RECREATE");

    float ecorebins[NECOREBIN + 1] = {0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 50.0};

    float  posbins[100];
    for (int i = 0; i < 100; i++)
        {
            posbins[i] = i;
        }

    for (Int_t ias = 0; ias < NARMSECT; ias++)
        {
            char hname[256];
            sprintf(hname, "hEmc3D_Armsect_%i", ias);

            int ny_max, nz_max;
            if (IsPbGl(ias))
                {
                    ny_max = YPOS_PBGL;
                    nz_max = ZPOS_PBGL;
                }
            else
                {
                    ny_max = YPOS_PBSC;
                    nz_max = ZPOS_PBSC;
                }

            hEmc3D[ias] = new TH3D(hname, "", nz_max, posbins, ny_max, posbins,
                                   NECOREBIN, ecorebins);
        }
    return 0;
}

int EMCalMapPP510::InitRun(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

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

int EMCalMapPP510::process_event(PHCompositeNode *topNode)
{
    // Informational message...
    nRunEvents++;

    if (nRunEvents % 1000 == 0 && verbosity)
        {
            if (verbo)
                {
                    cout << "Events for run " << runNumber << " = " << nRunEvents << endl;
                }
        }

    //  Get the data I need...
    PHGlobal *phglobal                     = getClass<PHGlobal>                (topNode, "PHGlobal");
    if (nRunEvents == 1 && !phglobal)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    float zvertex = phglobal->getBbcZVertex();
    if (fabs(zvertex) > 30.0)
        {
            return DISCARDEVENT;
        }

    emcClusterContainer *emcclustercontainer = getClass<emcClusterContainer>      (topNode, "emcClusterContainer");
    if (nRunEvents == 1 && !emcclustercontainer)
        {
            cout << "No emcClusterContainer!  No sense continuing" << endl;
            exit(1);
        }

    //The number of towers is 48(y) x 96(z) x 2(PbGl) + 36(y) x 72(z) x 6(PbSc) = 9216+15552 = 24768.
    //Tower ID is a unique ID of all towers in eight sectors (0~24767).

    int Nclus = emcclustercontainer->size();
    for (int iclus = 0; iclus < Nclus; iclus++)
        {
            emcClusterContent* clus = emcclustercontainer->getCluster(iclus);

            int arm  = clus->arm(); //In EMCal convention, West Arm is 0 and East Arm is 1, and thus armsector 0...7 are W0...W3 E0...E3
            int sector = clus->sector();
            int armsect = (arm * 4) + sector;
            int yPos = clus->iypos();
            int zPos = clus->izpos();
            float ecore = clus->ecore();

            hEmc3D[armsect]->Fill(zPos, yPos, ecore);
        }

    // any other return code might lead to aborting the event or analysis
    return 0;
}

int EMCalMapPP510::EndRun(const int runNumber)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    if (verbo)
        {
            cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
            cout << "Run number:                                 " << runNumber << endl;
            cout << "Number of total events for this run:        " << nRunEvents << endl;

        }
    return 0;
}

int EMCalMapPP510::End(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  End called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    outfile->Write(outfname.c_str());
    outfile->Close();

    if (verbo)
        {
            cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
        }
    return 0;
}


int EMCalMapPP510::IsPbGl(int armsect)
{
    return ((armsect == 4 || armsect == 5) ? 1 : 0);
}
