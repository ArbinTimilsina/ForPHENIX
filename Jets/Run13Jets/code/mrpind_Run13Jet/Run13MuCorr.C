#include <cmath>

#include <TOAD.h>
#include <TFile.h>
#include <TTree.h>
#include <Run13MuCorr.h>

using namespace std;

void ReadMuTree(int runnumber){

  mu = 0;

  cout << "LOADING MU FILE" << endl;

  TOAD* tofmaploader = new TOAD("Run13Jet");
  string locationtofmap = tofmaploader->location("muCorr/lumitree_run13Fin.root");
  TFile *tofmapf = new TFile(locationtofmap.c_str(), "read");

  if(!tofmapf->IsOpen())
    {
      cout << "Can not find " << locationtofmap.c_str() << "." << endl;
      delete tofmaploader;
      delete tofmapf;
      return;
    }
  TTree *Ttof = (TTree*)tofmapf->Get("T");
  if(Ttof==NULL){
    
    cout <<"Can not find T in tofmap root file." << endl;
    delete tofmaploader;
    delete tofmapf;
    
    return;
  }

  int runnum;
  double meannum;
  Ttof->SetBranchAddress("runnum", &runnum);
  Ttof->SetBranchAddress("mu", &meannum);

  int entries = Ttof->GetEntries();
  bool chk = false;
  for(int i = 0; i < entries; i++){//loop over tree

    Ttof->GetEntry(i);
    if(runnum == runnumber){

      mu = meannum;
      chk = true;
      break;

    }//endif
    
  }//end i

  if(!chk){
    cout << "ERROR: can't find Mu for this run" << endl;
    return;
  }

  delete tofmaploader;
  delete tofmapf;

  return;

}//end ReadMuTree()

double getMu(){

  double retval = 0;
  retval = mu;

  return retval;

}

