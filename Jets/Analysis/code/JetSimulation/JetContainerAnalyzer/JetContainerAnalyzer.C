//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <PHObject.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//PHPythia tools
#include <PHPythiaHeader.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TNtuple.h>

//For vertex during sim reco
#include <VtxOut.h>

//My source file
#include "JetContainerAnalyzer.h"

//For acceptance
#include "TrackQualityCuAu.h"
#include "TrackQualityPP.h"

using namespace std;
using namespace findNode;

typedef PHIODataNode<PHObject> PHObjectNode_t;

//================================ Constructor ================================
//Here we can initiate some variables
JetContainerAnalyzer::JetContainerAnalyzer(const float ir, const float iminPt,
					   const float incPythia, const float iminCfPythia, const float imaxCfPythia)
    : SubsysReco("JetContainerAnalyzer")
{
    R               = ir;
    minPt           = iminPt;
    ncPythia        = incPythia;
    minCfPythia     = iminCfPythia;
    maxCfPythia     = imaxCfPythia;

    nTotalEvents = 0;
    nGoodEvents = 0;

    eventContainer = NULL;
    jetContainer = NULL;
    particleContainer = NULL;

    return;
}

int JetContainerAnalyzer::Init(PHCompositeNode *topNode)
{

    //To access various nodes
    se = Fun4AllServer::instance();
    rc = recoConsts::instance();

    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");

    if(isCuAu)
        {
            jetAnalyzer->SetCuAu(true);
        }
    else
        {
            jetAnalyzer->SetCuAu(false);
        }
    jetAnalyzer->MyInit();

    //For PYTHIA trigger
    jetTriggerPythia = new JetTriggerPythia(R, ncPythia, minPt, minCfPythia, maxCfPythia, 1);

    CreateNodeTree(topNode);

    return EVENT_OK;
}

int JetContainerAnalyzer::InitRun(PHCompositeNode *topNode)
{

    return EVENT_OK;
}

int JetContainerAnalyzer::process_event(PHCompositeNode *topNode)
{
    nTotalEvents++;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Pythia info
    ///////////////////////////////////////////////////////////////////////////////////////////////
    pythiaNode = se->topNode(rc->get_CharFlag("PHPYTHIA_TOPNODE"));
    //Details about the event: event_gen/src/PHPythia/PHPythia/PHPythiaHeader.h
    PHPythiaHeader *phpythiaheader = findNode::getClass<PHPythiaHeader>(pythiaNode, "PHPythiaHeader");
    if (!phpythiaheader)
        {
            cout << "No PHPythiaHeader! No sense continuing" << endl;
            exit(1);
        }

    zVertexPythia = phpythiaheader->GetPrimaryVertexZ();
    zVertexPythia = zVertexPythia * MM2CM;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Pisa+Reco info
    ///////////////////////////////////////////////////////////////////////////////////////////////
    pisaRecoNode = se->topNode(rc->get_CharFlag("PISARECO_TOPNODE"));
    //For ZVertex used in reconstruction
    //https://www.phenix.bnl.gov/viewvc/viewvc.cgi/phenix/offline/packages/vtx/VtxOut.h?revision=1.25&view=co
    VtxOut *vertexOut = getClass<VtxOut>(pisaRecoNode, "VtxOut");
    if (!vertexOut)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    zVertexPisaReco = vertexOut->get_ZVertex();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Data info (only for Cu+Au)
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if(isCuAu)
        {
            dataNode = se->topNode(rc->get_CharFlag("DATA_TOPNODE"));

            phglobal = getClass<PHGlobal>(dataNode, "PHGlobal");
            if (!phglobal)
                {
                    cout << "No PHGlobal!  No sense continuing" << endl;
                    exit(1);
                }

            centrality = phglobal->getCentrality();
            zVertexData = phglobal->getBbcZVertex();

        }
    else
        {
            centrality = -9999.9;
            zVertexData = -9999.9;
        }

    eventContainer->SetBimpact(centrality);
    eventContainer->SetPrimaryVertexX(zVertexPythia);
    eventContainer->SetPrimaryVertexY(zVertexPisaReco);
    eventContainer->SetPrimaryVertexZ(zVertexData);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //Get True Jet
    ///////////////////////////////////////////////////////////////////////////////////////////////
    getPythiaTrueJet(pythia_particles, true_jets);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //After vertex feeding on the PYTHIA particles, true jet might not reconstruct
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if(true_jets.size() == 0)
        {
            return ABORTEVENT;
        }

    for(unsigned int t = 0; t < true_jets.size(); t++)
        {
	    TMCParticle tempJet;
	    tempJet.SetKF(true_jets[t].arm);
	    tempJet.SetPx(true_jets[t].pT);
	    tempJet.SetPy(true_jets[t].eta);
	    tempJet.SetPz(true_jets[t].phi);
	    tempJet.SetVx(true_jets[t].nc);
	    tempJet.SetVy(true_jets[t].cf);
	    tempJet.SetVz(true_jets[t].nf);
	    tempJet.SetTime(true_jets[t].disc);

	    jetContainer->addParticle(tempJet);
        }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    getPisaRecoParticles(pisa_reco_good_particles);
    getTotalParticles(pisa_reco_good_particles, total_particles);

    for(unsigned int p = 0; p < total_particles.size(); p++)
        {
	    TMCParticle tempParticle;
	    tempParticle.SetKF(total_particles[p].arm);
	    tempParticle.SetKS(total_particles[p].charge);
	    tempParticle.SetEnergy(total_particles[p].energy);
	    tempParticle.SetVx(total_particles[p].mom);
	    tempParticle.SetVy(total_particles[p].pT);
	    tempParticle.SetVz(total_particles[p].eT);
	    tempParticle.SetPx(total_particles[p].px);
	    tempParticle.SetPy(total_particles[p].py);
	    tempParticle.SetPz(total_particles[p].pz);
	    tempParticle.SetMass(total_particles[p].eta);
	    tempParticle.SetTime(total_particles[p].phi);
	    tempParticle.SetLifetime(total_particles[p].phiDC);

	    particleContainer->addParticle(tempParticle);
        }

    nGoodEvents++;

    if(nGoodEvents == 10)
        {
            cout << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Vertex info for 10th event: " << endl;
            cout << "PYTHIA: " << zVertexPythia << endl;
            cout << "PISA + RECO: " << zVertexPisaReco << endl;
            cout << "Data: " << zVertexData << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << endl;
        }

    return 0;
}

void JetContainerAnalyzer::SetCuAu(bool what = true)
{
    isCuAu = what;
}

void JetContainerAnalyzer::getPythiaTrueJet(std::vector<particles>& particle_list,
					    std::vector<jets>& true_jet_list)
{
    jetTriggerPythia->GetPythia(pythiaNode, particle_list, true_jet_list);
}


void JetContainerAnalyzer::getPisaRecoParticles(std::vector<particles>& good_particle_list)
{
    //Get Tracks
    std::vector<tracks> all_track_list;
    all_track_list.clear();
    jetAnalyzer->GetTracks(pisaRecoNode, all_track_list);

    //Get Clusters
    std::vector<clusters> all_cluster_list;
    all_cluster_list.clear();
    jetAnalyzer->SetData(false);
    jetAnalyzer->GetClusters(pisaRecoNode, all_cluster_list, zVertexPisaReco);

    std::vector<tracks> good_track_list;
    good_track_list.clear();
    std::vector<clusters> good_cluster_list;
    good_cluster_list.clear();
    jetAnalyzer->GetParticles(all_track_list, all_cluster_list, good_track_list, good_cluster_list, good_particle_list);
}

void JetContainerAnalyzer::getTotalParticles(std::vector<particles> sim_good_particle_list,
					     std::vector<particles>& total_particle_list)
{
    //Get Tracks
    std::vector<tracks> cnt_track_list;
    cnt_track_list.clear();
    if(isCuAu)
        {
            jetAnalyzer->GetTracks(dataNode, cnt_track_list);
        }

    //Get Clusters
    std::vector<clusters> cnt_cluster_list;
    cnt_cluster_list.clear();
    if(isCuAu)
        {
            jetAnalyzer->SetData(true);
            jetAnalyzer->GetClusters(dataNode, cnt_cluster_list, zVertexData);
        }

    std::vector<tracks> cnt_good_track_list;
    cnt_good_track_list.clear();
    std::vector<clusters> cnt_good_cluster_list;
    cnt_good_cluster_list.clear();
    std::vector<particles> cnt_good_particle_list;
    cnt_good_particle_list.clear();
    if(isCuAu)
        {
            jetAnalyzer->GetParticles(cnt_track_list, cnt_cluster_list, cnt_good_track_list, cnt_good_cluster_list, cnt_good_particle_list);
        }

    //Do the embedding- tracks and clusters
    total_particle_list.reserve(sim_good_particle_list.size() + cnt_good_particle_list.size());
    total_particle_list.insert(total_particle_list.end(), sim_good_particle_list.begin(), sim_good_particle_list.end());
    total_particle_list.insert(total_particle_list.end(), cnt_good_particle_list.begin(), cnt_good_particle_list.end());
}


int JetContainerAnalyzer::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetContainerAnalyzer::End(PHCompositeNode * topNode)
{
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nTotalEvents << endl;
    cout << "Total events with true jet = " << nGoodEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    return EVENT_OK;
}

int JetContainerAnalyzer::CreateNodeTree(PHCompositeNode * topNode)
{
    //*******************************************************************************************************
    //Manage nodes
    //*******************************************************************************************************
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
        {
            cout << PHWHERE << "DST Node missing doing nothing" << endl;
            return ABORTRUN;
        }

    //Event information
    eventContainer = new PHHijingHeaderV3();
    PHObjectNode_t *EventContainer = new PHObjectNode_t(eventContainer, "EventContainerNode", "PHObject");
    dstNode->addNode(EventContainer);

    //True Jet information
    jetContainer = new PHPythiaContainerV2();
    PHObjectNode_t *JetContainer = new PHObjectNode_t(jetContainer, "JetContainerNode", "PHObject");
    dstNode->addNode(JetContainer);

    //Tracks/Clusters information
    particleContainer = new PHPythiaContainerV2();
    PHObjectNode_t *ParticleContainer = new PHObjectNode_t(particleContainer, "ParticleContainerNode", "PHObject");
    dstNode->addNode(ParticleContainer);

    return 0;
}

int JetContainerAnalyzer::ResetEvent(PHCompositeNode *topNode)
{
    pythia_particles.clear();
    true_jets.clear();

    pisa_reco_good_particles.clear();
    total_particles.clear();

    PHNodeIterator mainIter(topNode);
    PHNodeReset reset;
    if (mainIter.cd(Name()))
        {
            mainIter.forEach(reset);
        }

    return 0;
}
























































