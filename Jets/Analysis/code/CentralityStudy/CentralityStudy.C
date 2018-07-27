//General PHENIX tools
#include <getClass.h>
#include <PHCompositeNode.h>
#include <phool.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TNtuple.h>

#include <VtxOut.h>

#include <stdlib.h>  

//My source file
#include "CentralityStudy.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
CentralityStudy::CentralityStudy()
  : SubsysReco("CentralityStudy")
{
  nTotalEvents = 0;

  return;
}

int CentralityStudy::Init(PHCompositeNode *topNode)
{
  //Output file name
  outfile = new TFile("CentralityStudy.root", "RECREATE");

  //Histograms
  hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);
  hVertex = new TH1F("hVertex", "Vertex", 100, -50, 50);
  hBbcQtotal = new TH1F("hBbcQtotal","BBC q total", 10000, 0.0, 2500.0);

  //*******************************************************************************************************
  return EVENT_OK;
}

int CentralityStudy::InitRun(PHCompositeNode *topNode)
{
  return EVENT_OK;
}

int CentralityStudy::ResetEvent(PHCompositeNode *topNode)
{
  return EVENT_OK;
}


int CentralityStudy::process_event(PHCompositeNode *topNode)
{
  nTotalEvents++;

  VtxOut *vertexOut = getClass<VtxOut>(topNode, "VtxOut");
  if (!vertexOut)
      {
	  cout << "No VtxOut!  No sense continuing" << endl;
	  exit(1);
      }

  float zVertexOut = vertexOut->get_ZVertex();

  if(nTotalEvents ==1){
      cout<<"Vertex from VtxOut is: "<<zVertexOut<<" cm!!!"<<endl;
      }


  PHGlobal *phGlobal = getClass<PHGlobal>(topNode, "PHGlobal");
  if (!phGlobal)
      {
          cout << "No PHGlobal!  No sense continuing" << endl;
          exit(1);
      }

  float zVertexBBC  = phGlobal->getBbcZVertex();
  if(nTotalEvents ==1){
      cout<<"Vertex from PHGlobal is: "<<zVertexBBC<<" cm!!!"<<endl<<endl;
  }
  hVertex->Fill(zVertexBBC);

  float BbcQs = phGlobal->getBbcChargeS();
  float BbcQn = phGlobal->getBbcChargeN();
  float BbcQtotal = BbcQs + BbcQn;
  hBbcQtotal->Fill(BbcQtotal);

  //int BbcNs = phGlobal->getBbcMultS();
  //int BbcNn = phGlobal->getBbcMultN();
  //int BbcNtotal = BbcNs + BbcNn;

  //float centrality = phGlobal->getCentrality();

  return 0;
}



int CentralityStudy::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int CentralityStudy::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nTotalEvents);

    cout << endl;
    cout << "Total events processed: " << nTotalEvents << endl;
    outfile->Write();
    outfile->Close();

    return EVENT_OK;
}

