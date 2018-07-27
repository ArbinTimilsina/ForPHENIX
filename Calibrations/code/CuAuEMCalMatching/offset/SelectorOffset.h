#ifndef SelectorOffset_h
#define SelectorOffset_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TMath.h>
#include <TH2.h>

static const int NSECT = 8;
static const int NPT = 14;
static const int NTHETA = 8;

class SelectorOffset : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Float_t         armsect;
   Float_t         zed;
   Float_t         centrality;
   Float_t         pT;
   Float_t         charge;
   Float_t         phi;
   Float_t         emcdphi;
   Float_t         emcdz;
   Float_t         theta;
   Float_t         beta;

   // List of branches
   TBranch        *b_armsect;   //!
   TBranch        *b_zed;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_pT;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_emcdphi;   //!
   TBranch        *b_emcdz;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_beta;   //!

   std::string rootfname;
   unsigned long long  eventcounter;

   TH1F *hEmcdZ[NSECT][NPT][NTHETA];
   TH2F *hEmcBefore[NSECT][NPT];
   TH2F *hEmcAfter[NSECT][NPT];

   float slope[NSECT][NPT];

   int GetPtBin(float);
   void LoadSlope(const char* filename);

   SelectorOffset(TTree * /*tree*/ =0) { }
   virtual ~SelectorOffset() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(SelectorOffset,0);
};

#endif

#ifdef SelectorOffset_cxx
void SelectorOffset::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("armsect", &armsect, &b_armsect);
   fChain->SetBranchAddress("zed", &zed, &b_zed);
   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("pT", &pT, &b_pT);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("emcdphi", &emcdphi, &b_emcdphi);
   fChain->SetBranchAddress("emcdz", &emcdz, &b_emcdz);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
}

Bool_t SelectorOffset::Notify()
{
   return kTRUE;
}

#endif // #ifdef SelectorOffset_cxx
