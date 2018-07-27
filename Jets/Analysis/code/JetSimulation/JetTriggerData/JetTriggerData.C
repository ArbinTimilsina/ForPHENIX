//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//Root tools
#include <TFile.h>
#include <TNtuple.h>

//My source file
#include "JetTriggerData.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
JetTriggerData::JetTriggerData(const float ir, const float inc, const float iminPt, const float iminCf, const float imaxCf,
                               const int iTrigEvents)
    : SubsysReco("JetTriggerData")
{
    R           = ir;
    nc          = inc;
    minPt       = iminPt;
    minCf       = iminCf;
    maxCf       = imaxCf;

    nTrigEvents = iTrigEvents;

    nTotalEvents                  = 0;
    nTotalFailedVertexEvents      = 0;
    nTotalJetEvents               = 0;
    nTriggeredEvents              = 0;

    return;
}

int JetTriggerData::Init(PHCompositeNode *topNode)
{
    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(true);
    jetAnalyzer->SetCuAu(true);
    jetAnalyzer->MyInit();

    vertexFile.open("vertex.txt");

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Triggering on events without jet with following conditon:" << endl;
    cout << "R = " << R << endl;
    cout << "nc >= " << nc << endl;
    cout << "min pT = " << minPt << " (GeV/c)" << endl;
    cout << "cf between: " << minCf << " to " << maxCf << endl;
    cout << "no of triggered events required = " << nTrigEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;


    return EVENT_OK;
}

int JetTriggerData::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int JetTriggerData::ResetEvent(PHCompositeNode *topNode)
{
    return EVENT_OK;
}


int JetTriggerData::process_event(PHCompositeNode *topNode)
{
    if(nTriggeredEvents >= nTrigEvents)
        {
            return ABORTRUN;
        }


    nTotalEvents++;

    PHGlobal *phglobal = getClass<PHGlobal>(topNode, "PHGlobal");
    if (!phglobal)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    EventHeader *evtheader = getClass<EventHeader>(topNode, "EventHeader");
    if (!evtheader)
        {
            cout << "No EventHeader! No sense continuing" << endl;
            exit(1);
        }

    int eventNumber = evtheader->get_EvtSequence();

    if(nTotalEvents == 1)
        {
            cout << endl;
            cout << "////////////////////////////////////////////" << endl;
            cout << "First event number is: " << eventNumber << endl;
            cout << "////////////////////////////////////////////" << endl << endl;
        }

    zvertex = phglobal->getBbcZVertex();
    centrality = phglobal->getCentrality();

    bool failVertex = (fabs(zvertex) > VERTEX_CUT);

    if (failVertex)
        {
            nTotalFailedVertexEvents++;
            return ABORTEVENT;
        }

    std::vector<tracks> all_tracks;
    all_tracks.clear();
    std::vector<clusters> all_clusters;
    all_clusters.clear();
    std::vector<tracks> charged_particles;
    charged_particles.clear();
    std::vector<clusters> neutral_particles;
    neutral_particles.clear();
    std::vector<particles> all_particles;
    all_particles.clear();

    std::vector<jets> antikt_jets;
    antikt_jets.clear();
    std::vector<particles> jet_constituents;
    jet_constituents.clear();

    //Get Tracks
    jetAnalyzer->GetTracks(topNode, all_tracks);

    //Get Clusters
    jetAnalyzer->GetClusters(topNode, all_clusters, zvertex);

    jetAnalyzer->GetParticles(all_tracks, all_clusters, charged_particles, neutral_particles, all_particles);

    jetAnalyzer->GetAntiKtCommon(all_particles, R, nc, minPt, minCf, maxCf, antikt_jets, jet_constituents,
                                 -1);

    bool triggered = false;
    if(antikt_jets.size() == 0)
        {
            triggered = true;
        }

    if(triggered)
        {
            //Write out vertex file
            vertexFile << nTriggeredEvents << " " << zvertex << " " << centrality << "\n";

            nTriggeredEvents++;
            return EVENT_OK;
        }
    else
        {
            nTotalJetEvents++;
            return ABORTEVENT;
        }
}

int JetTriggerData::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetTriggerData::End(PHCompositeNode *topNode)
{
    cout << endl << endl;
    cout << "Total events processed:                  " << nTotalEvents << endl;
    cout << "Total events that failed vertex cut:     " << nTotalFailedVertexEvents << endl;
    cout << "Total events that had jet:               " << nTotalJetEvents << endl;
    cout << "Total triggered events:                  " << nTriggeredEvents << endl << endl;

    vertexFile.close();

    return EVENT_OK;
}






















