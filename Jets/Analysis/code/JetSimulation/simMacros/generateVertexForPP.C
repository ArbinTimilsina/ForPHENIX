#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <TNtuple.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TChain.h>
#include <sstream>
#include <TRandom3.h>

using namespace std;

void generateVertexForPP(int nevents = 10)
{
  ofstream vertexFile;
  vertexFile.open("vertex.txt");

  TH1F *hVertex = new TH1F("hVertex", "Gaussian p+p MB vertex distribution", 500, -50, 50);

  int events = 0;
  while(events < nevents)
    {
      TRandom3 *myRandom3 = new TRandom3();
      myRandom3->SetSeed(0);
      float zvertex = myRandom3->Gaus(0.0, 15);
      hVertex->Fill(zvertex);

      if(fabs(zvertex) <= 10.0)
        {
	  vertexFile << events << " " << zvertex << " " << 0.0 << "\n";
	  events++;
        }
    }
  vertexFile.close();

  TCanvas *c = new TCanvas;
  hVertex->Draw();
  c->Update();
  c->SaveAs("PythiaVertex.png");
  c->Close();
}


