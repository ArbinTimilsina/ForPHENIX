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
#include <PHPyTrigger.h>
#include <PHPythiaHeader.h>
#include <PHPyCommon.h>
#include <PHPythiaContainer.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

//Root tools
#include <TH1.h>
#include <TFile.h>
#include <TNtuple.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//My source file
#include "Run13JetTriggerPythia.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
Run13JetTriggerPythia::Run13JetTriggerPythia(const float iminPt, const int iTrigEvents)
    : SubsysReco("Run13JetTriggerPythia")
{
    minPt        = iminPt;
    nTrigEvents  = iTrigEvents;

    nTotalEvents      = 0;
    nTriggeredEvents  = 0;
    nTotalJets        = 0;

    nTotalParticles   = 0;

    return;
}

int Run13JetTriggerPythia::Init(PHCompositeNode *topNode)
{
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Triggering on the following conditions for PYTHIA jet:" << endl;
    cout << "min pT = " << minPt << " (GeV/c)" << endl;
    cout << "no of triggered events required = " << nTrigEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    return EVENT_OK;
}

int Run13JetTriggerPythia::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int Run13JetTriggerPythia::ResetEvent(PHCompositeNode *topNode)
{
    pythia_particles.clear();
    true_jets.clear();

    return EVENT_OK;
}


int Run13JetTriggerPythia::process_event(PHCompositeNode *topNode)
{
    if(nTriggeredEvents >= nTrigEvents)
        {
            return ABORTRUN;
        }

    //Get Pythia
    GetPythia(topNode, pythia_particles, true_jets);

    nTotalEvents++;

    bool triggered = false;
    if(true_jets.size() > 0)
        {
            triggered = true;
        }

    if(triggered)
        {
            nTriggeredEvents++;
            nTotalJets += true_jets.size();
            nTotalParticles += pythia_particles.size();

            return EVENT_OK;
        }
    else
        {
            return ABORTEVENT;
        }
}

void Run13JetTriggerPythia::GetPythia(PHCompositeNode *topNode,
                                      std::vector<pythiaParticles>& particle_list,
                                      std::vector<pythiaJets>& jet_list)
{
    PHPythiaContainer *phpythia = getClass<PHPythiaContainer> (topNode, "PHPythia");
    if (!phpythia)
        {
            cout << "No PHPythia! No sense continuing" << endl;
            exit(1);
        }

    int npart = phpythia->size();
    for (int ipart = 0; ipart < npart; ipart++)
        {
            TMCParticle *part = phpythia->getParticle(ipart);

            //Look for final state particles
            bool passKS = (part->GetKS() == 1);

            //Ignore all beam remnants
            bool passBeanRemnants = (sqrt(pow(part->GetPx(), 2) + pow(part->GetPy(), 2)) != 0.0);

            if(passKS && passBeanRemnants)
                {
                    //KF- falvor
                    int kf = part->GetKF();

                    float energy = part->GetEnergy();
                    TLorentzVector *fVector = new TLorentzVector(part->GetPx(), part->GetPy(), part->GetPz(), energy);
                    float eta   = fVector->Eta();
                    float theta = fVector->Theta();
                    float mom   = fVector->P();
                    float pT    = mom * sin(theta);
                    float phi   = phiReduce(fVector->Phi());

                    TPythia6 *tpythia6 = new TPythia6();
                    int charge = tpythia6->Pychge(kf);
                    if (charge == 3)
                        {
                            charge = 1;
                        }
                    if (charge == -3)
                        {
                            charge = -1;
                        }
                    delete tpythia6;

                    pythiaParticles temp;
                    temp.charge      = charge;

                    temp.energy      = energy;
                    temp.mom         = mom;
                    temp.pT          = pT;
                    temp.px          = part->GetPx();
                    temp.py          = part->GetPy();
                    temp.pz          = part->GetPz();
                    temp.eta         = eta;
                    temp.phi         = phi;

                    //No neutrino
                    bool passNoNeutrino = ((kf != PY_NU_E) && (kf != -PY_NU_E) &&
                                           (kf != PY_NU_MU) && (kf != -PY_NU_MU) &&
                                           (kf != PY_NU_TAU) && (kf != -PY_NU_TAU));

                    //No muons
                    bool passNoMuon = ((kf != PY_MU) && (kf != -PY_MU));
                    bool passEverything = passNoNeutrino && passNoMuon;

                    if(passEverything)
                        {
                            particle_list.push_back(temp);
                        }
                }
        }

    //Make jet
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.3, fastjet::E_scheme, fastjet::Best);

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

            pythiaJets tempJet;
            tempJet.pT            = aFastJet.perp();
            tempJet.eta           = aFastJet.pseudorapidity();
            tempJet.phi           = phiReduce(aFastJet.phi());

            if((fabs(tempJet.eta) < 0.35) && (tempJet.pT > minPt))
                {
                    jet_list.push_back(tempJet);
                }
        }
    delete antikt;

}

int Run13JetTriggerPythia::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int Run13JetTriggerPythia::End(PHCompositeNode *topNode)
{
    cout << endl;
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed:               " << nTotalEvents << endl;
    cout << "Total triggered events:               " << nTriggeredEvents << endl;
    cout << "Total Jets found:                     " << nTotalJets << endl;
    cout << "Total particles:                      " << nTotalParticles << endl;
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    return EVENT_OK;
}
