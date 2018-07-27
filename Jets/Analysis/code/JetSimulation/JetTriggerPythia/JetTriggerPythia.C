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

//My source file
#include "JetTriggerPythia.h"

using namespace std;
using namespace findNode;

//================================ Constructor ================================
//Here we can initiate some variables
JetTriggerPythia::JetTriggerPythia(const float ir, const float inc, const float iminPt, const float iminCf, const float imaxCf,
                                   const int iTrigEvents)
    : SubsysReco("JetTriggerPythia")
{
    R            = ir;
    nc           = inc;
    minPt        = iminPt;
    minCf        = iminCf;
    maxCf        = imaxCf;

    nTrigEvents  = iTrigEvents;

    nTotalEvents      = 0;
    nTriggeredEvents  = 0;
    nTotalJets        = 0;

    nTotalParticles   = 0;

    return;
}

int JetTriggerPythia::Init(PHCompositeNode *topNode)
{
    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(false);

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Triggering on the following conditions for PYTHIA jet:" << endl;
    cout << "R = " << R << endl;
    cout << "nc >= " << nc << endl;
    cout << "min pT = " << minPt << " (GeV/c)" << endl;
    cout << "cf between: " << minCf << " to " << maxCf << endl;
    cout << "no of triggered events required = " << nTrigEvents << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    return EVENT_OK;
}

int JetTriggerPythia::InitRun(PHCompositeNode *topNode)
{
    return EVENT_OK;
}

int JetTriggerPythia::ResetEvent(PHCompositeNode *topNode)
{
    pythia_particles.clear();
    true_jets.clear();

    return EVENT_OK;
}


int JetTriggerPythia::process_event(PHCompositeNode *topNode)
{
   if(nTriggeredEvents >= nTrigEvents)
        {
            return ABORTRUN;
        }

    GetPythia(topNode, pythia_particles, true_jets);

    nTotalEvents++;

    bool triggered = false;
    if(true_jets.size() > 0)
        {
            triggered = true;
        }

    //cout<<"Total events: "<<nTotalEvents<<endl;
    //cout<<"Triggered events: "<<nTriggeredEvents<<endl;
    //cout<<"Jet size: "<<true_jets.size()<<endl;
    //cout<<"Particle size: "<<pythia_particles.size()<<endl<<endl;
    if(triggered)
        {
            nTriggeredEvents++;
            nTotalJets += true_jets.size();
            nTotalParticles += pythia_particles.size();

            if(nTriggeredEvents == 1)
                {
                    cout << endl;
                    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                    cout << "For first triggered event: " << endl;
                    cout << "Particle's:" << setw(15) << "charge" << setw(15) << "energy" << setw(15) << "mom" << setw(15)
                         << "pT" << setw(15) << "px" << setw(15) << "py" << setw(15) << "pz" << setw(15)
                         << "eta" << setw(15) << "phi" << endl;
                    for(int p = 0; p < pythia_particles.size(); p++)
                        {
                            cout << "Particle #: " << p + 1 << setw(15) << pythia_particles[p].charge << setw(15)
                                 << pythia_particles[p].energy << setw(15) << pythia_particles[p].mom << setw(15)
                                 << pythia_particles[p].pT << setw(15) << pythia_particles[p].px << setw(15)
                                 << pythia_particles[p].py << setw(15) << pythia_particles[p].pz << setw(15)
                                 << pythia_particles[p].eta << setw(15) << pythia_particles[p].phi << endl;
                        }
                    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

                    cout << endl;
                    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                    cout << "True Jet Reconstructed with: " << endl;
                    cout << "nc: " << true_jets[0].nc << ", pT: " << true_jets[0].pT << ", cf: " << true_jets[0].cf << ", eta: " <<
			true_jets[0].eta << ", phi: " << true_jets[0].phi << endl;
                    cout << "***********************************************************************" << endl << endl;
                }
            return EVENT_OK;
        }
    else
        {
            return ABORTEVENT;
        }
}


void JetTriggerPythia::GetPythia(PHCompositeNode *topNode,
                                 std::vector<particles>& particle_list,
                                 std::vector<jets>& antikt_jet_list)
{
    PHPythiaContainer *phpythia = getClass<PHPythiaContainer> (topNode, "PHPythia");
    if (!phpythia)
        {
            cout << "No PHPythia! No sense continuing" << endl;
            exit(1);
        }

    PHPythiaHeader *phpythiaheader = findNode::getClass<PHPythiaHeader>(topNode, "PHPythiaHeader");
    if (!phpythiaheader)
        {
            cout << "No PHPythiaHeader! No sense continuing" << endl;
            exit(1);
        }

    double zvertex = phpythiaheader->GetPrimaryVertexZ();
    zvertex = zvertex * 0.1;

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
                    float eT    = energy * sin(theta);
                    float phi   = jetAnalyzer->phiReduce(fVector->Phi());

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
                    temp.phiDC       = -999.9;

                    temp.ertTrigger  = false;

                    temp.id          = kf;
                    temp.jetPt       = -999.9;

                    //No neutrino
                    bool passNoNeutrino = ((kf != PY_NU_E) && (kf != -PY_NU_E) &&
                                           (kf != PY_NU_MU) && (kf != -PY_NU_MU) &&
                                           (kf != PY_NU_TAU) && (kf != -PY_NU_TAU));

                    //No muons
                    bool passNoMuon = ((kf != PY_MU) && (kf != -PY_MU));

                    //Particle should pass all the criteria
                    bool passEverything = passNoNeutrino && passNoMuon;

                    if(passEverything)
                        {
                            particle_list.push_back(temp);
                        }
                }
        }
    std::vector<jets> true_jet_list;
    true_jet_list.clear();
    std::vector<particles> constituent_list;
    constituent_list.clear();
    jetAnalyzer->GetAntiKtCommon(particle_list, R, nc, minPt, minCf, maxCf, true_jet_list, constituent_list, -1);

    //Make only eta cut
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
	    float jetEta = true_jet_list[t].eta;
	    bool satisfyEta = fabs(jetEta) < 0.35;
            if(satisfyEta)
                {
                    antikt_jet_list.push_back(true_jet_list[t]);
                }
        }
}

int JetTriggerPythia::EndRun(const int runNumber)
{
    return EVENT_OK;
}

int JetTriggerPythia::End(PHCompositeNode *topNode)
{
    cout << endl;
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed:               " << nTotalEvents << endl;
    cout << "Total triggered events:               " << nTriggeredEvents << endl;
    cout << "Total Jets found:                     " << nTotalJets << endl;
    cout << "Total particles in PHENIX acceptance: " << nTotalParticles << endl;
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    return EVENT_OK;
}






































