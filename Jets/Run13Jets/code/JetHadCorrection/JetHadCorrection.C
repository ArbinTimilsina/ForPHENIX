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

//PHPythia tools
#include <PHPythiaHeader.h>
#include <PHPyCommon.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//C tools
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <algorithm> // remove and remove_if

//Root tools
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TNtuple.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//My source file
#include "JetHadCorrection.h"

using namespace std;
using namespace findNode;
using namespace fastjet;

//================================ Constructor ================================
//Here we can initiate some variables
JetHadCorrection::JetHadCorrection()
    : SubsysReco("JetHadCorrection")
{
    writeTree = false;

    return;
}

int JetHadCorrection::Init(PHCompositeNode *topNode)
{
    //Output file name
    outfile = new TFile("JetHadCorrection.root", "RECREATE");
    nRunEvents = 0;

    //Histograms
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    hPtPartonJetsR2 = new TH1F("hPtPartonJetsR2", "p_{T} of Parton Jets, R = 0.2", NPTBINS, PTBINS);
    hPtParticleJetsR2 = new TH1F("hPtParticleJetsR2", "p_{T} of Particle Jets, R = 0.2", NPTBINS, PTBINS);

    hPtPartonJetsR3 = new TH1F("hPtPartonJetsR3", "p_{T} of Parton Jets, R = 0.3", NPTBINS, PTBINS);
    hPtParticleJetsR3 = new TH1F("hPtParticleJetsR3", "p_{T} of Particle Jets, R = 0.3", NPTBINS, PTBINS);

    //nTuples
    if(writeTree)
        {
            pythiaParton = new TTree("pythiaParton", "Partons of Pythia");
            pythiaParton->Branch("pythiaParton", &pythia_parton_list);

            pythiaParticle = new TTree("pythiaParticle", "Particles of Pythia");
            pythiaParticle->Branch("pythiaParticle", &pythia_particle_list);
        }

    return 0;
}

int JetHadCorrection::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int JetHadCorrection::ResetEvent(PHCompositeNode *topNode)
{
    pythia_parton_list.clear();
    pythia_particle_list.clear();

    return EVENT_OK;
}


int JetHadCorrection::process_event(PHCompositeNode *topNode)
{
    nRunEvents++;
    if (nRunEvents % 1000 == 0)
        {
            cout << "Events = " << nRunEvents << endl;
        }

    //*****************************************************************************************************************************
    //Get PYTHIA stuff
    //*****************************************************************************************************************************
    phpythia = getClass<PHPythiaContainer> (topNode, "PHPythia");
    if (!phpythia)
        {
            cout << "No PHPythia! No sense continuing" << endl;
            exit(1);
        }

    int npart = phpythia->size();
    for (int ipart = 0; ipart < npart; ipart++)
        {
            TMCParticle *part = phpythia->getParticle(ipart);

            //KF- falvor
            int kf = part->GetKF();

            float energy = part->GetEnergy();
            TLorentzVector *fVector = new TLorentzVector(part->GetPx(), part->GetPy(), part->GetPz(), energy);
            float eta   = fVector->Eta();
            float theta = fVector->Theta();
            float mom   = fVector->P();
            float pT    = mom * sin(theta);
            float eT    = energy * sin(theta);
            float phi   = phiReduce(fVector->Phi());

            TPythia6 *tpythia6 = new TPythia6();
            int charge = tpythia6->Pychge(kf);
            delete tpythia6;

            int arm = 1;
            if (phi > 1.57)
                {
                    arm = 0;
                }

            particles temp;
            temp.arm         = arm;
            temp.charge      = charge;
            temp.energy      = energy;
            temp.mom         = mom;
            temp.pT          = pT;
            temp.eT          = eT;
            temp.px          = part->GetPx();
            temp.py          = part->GetPy();
            temp.pz          = part->GetPz();
            temp.eta         = eta;
            temp.phi         = phi;

            //For partons, save only collision products
            //Note: For these, MPI and hadronization are turned off in pythia.cfg
            bool passParton = ((phpythia->getLineNumber(part) > 8) && (part->GetParent() <= 8));
            if(passParton)
                {
                    pythia_parton_list.push_back(temp);
                }

            //For particles, save final state particles
            bool passParticle = (part->GetKS() == 1);
            if(passParticle)
                {
                    if (charge == 3)
                        {
                            charge = 1;
                        }
                    if (charge == -3)
                        {
                            charge = -1;
                        }
                    temp.charge = charge;
                    pythia_particle_list.push_back(temp);
                }
        }

    //*****************************************************************************************************************************
    //Make jets
    //*****************************************************************************************************************************
    GetAntiKtJets(pythia_parton_list, 0.2, hPtPartonJetsR2);
    GetAntiKtJets(pythia_particle_list, 0.2, hPtParticleJetsR2);

    GetAntiKtJets(pythia_parton_list, 0.3, hPtPartonJetsR3);
    GetAntiKtJets(pythia_particle_list, 0.3, hPtParticleJetsR3);

    if(writeTree)
        {
            pythiaParton->Fill();
            pythiaParticle->Fill();
        }

    // any other return code might lead to aborting the event or analysis
    return 0;
}

void JetHadCorrection::GetAntiKtJets(std::vector<particles> particle_list, float R, TH1F *hJet)
{
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    //Generate input for jet reconstuction
    std::vector<fastjet::PseudoJet> jetParticles_all;
    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            fastjet::PseudoJet pseudoParticle(particle_list[h].px,
                                              particle_list[h].py,
                                              particle_list[h].pz,
                                              particle_list[h].energy);

            jetParticles_all.push_back(pseudoParticle);
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float jetEta = aFastJet.pseudorapidity();
            if(fabs(jetEta) < 0.35)
                {
                    hJet->Fill(aFastJet.perp());
                }
        }
    delete antikt;
}

int JetHadCorrection::End(PHCompositeNode *topNode)
{
    hEvents->SetBinContent(1, nRunEvents);

    outfile->Write();
    outfile->Close();

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nRunEvents << endl;

    return 0;
}


















