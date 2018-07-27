//General PHENIX tools
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHIODataNode.h>
#include <phool.h>
#include <getClass.h>
#include <RunHeader.h>

//HepMc tools
#include <HepMC/GenEvent.h>

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
#include "JetNloCorrection.h"

using namespace std;
using namespace findNode;
using namespace fastjet;

//================================ Constructor ================================
//Here we can initiate some variables
JetNloCorrection::JetNloCorrection()
    : SubsysReco("JetNloCorrection")
{
    return;
}

int JetNloCorrection::Init(PHCompositeNode *topNode)
{
    //For data code
    jetAnalyzer = new JetAnalyzer("dummyFile.root");
    jetAnalyzer->SetData(false);

    //Output file name
    outfile = new TFile("JetNloCorrection.root", "RECREATE");

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //Histograms
    //*****************************************************************************************************************************
    //General
    hEvents = new TH1F("hEvents", "Number of events", 2, 0, 2);

    hPtTrueJet = new TH1F("hPtTrueJet", "p_{T} of True Jets",
                          NPTBINS_FINAL, PTBINS_FINAL);

    hPtRecoJet = new TH1F("hPtRecoJet", "p_{T} of Reco Jets",
                          NPTBINS_FINAL, PTBINS_FINAL);

    hDistance = new TH1F("hDistance", "#DeltaR between True Jets and Reco Jets",
                         400, 0.0, 0.4);

    hPtMatchedJet = new TH1F("hPtMatchedJet", "p_{T} of Matched Jets",
                             NPTBINS_FINAL, PTBINS_FINAL);

    hPtForEfficiency = new TH1F("hPtForEfficiency", "p_{T} for Matching Efficiency",
                                NPTBINS_FINAL, PTBINS_FINAL);

    hResponseMatrix = new TH2F("hResponseMatrix", "p_{T, matched} vs p_{T, true} for Reco Jets",
                               NPTBINS_FINAL, PTBINS_FINAL, NPTBINS_FINAL, PTBINS_FINAL);


    for(unsigned int c = 0; c < 17; c++)
        {
            hForCorrectionFactor[c] = new TH1F(Form("hForCorrectionFactor_%u", c), "For Correction Factor",
                                               60.0, 0.0, 2.0);
        }

    //*****************************************************************************************************************************
    //*****************************************************************************************************************************
    //nTuples
    //*****************************************************************************************************************************
    pythiaParton = new TTree("pythiaParton", "Partons of Pythia");
    pythiaParton->Branch("pythiaParton", &pythia_parton_list);

    pythiaFinalParticle = new TTree("pythiaFinalParticle", "Final particles of Pythia");
    pythiaFinalParticle->Branch("pythiaFinalParticle", &pythia_final_particle_list);


    pythiaPartonJet = new TTree("pythiaPartonJet", "Jet reconstructed from partons of Pythia");
    pythiaPartonJet->Branch("pythiaPartonJet", &pythia_parton_jet_list);

    pythiaFinalParticleJet = new TTree("pythiaFinalParticleJet", "Jet reconstructed form final particles of Pythia");
    pythiaFinalParticleJet->Branch("pythiaFinalParticleJet", &pythia_final_particle_jet_list);

    pythiaMatchedJet = new TTree("pythiaMatchedJet", "Matched jet of Pythia");
    pythiaMatchedJet->Branch("pythiaMatchedJet", &pythia_matched_jet_list);
    //*****************************************************************************************************************************
    return 0;
}

int JetNloCorrection::InitRun(PHCompositeNode *topNode)
{
    nRunEvents = 0;

    return EVENT_OK;
}

int JetNloCorrection::ResetEvent(PHCompositeNode *topNode)
{
    pythia_parton_list.clear();
    pythia_final_particle_list.clear();

    pythia_parton_jet_list.clear();
    pythia_final_particle_jet_list.clear();

    pythia_matched_jet_list.clear();

    return EVENT_OK;
}


int JetNloCorrection::process_event(PHCompositeNode *topNode)
{
    // Informational message...
    nRunEvents++;
    if (nRunEvents % 1000 == 0)
        {
            cout << "Events = " << nRunEvents << endl;
        }

    //*****************************************************************************************************************************
    //Pythia
    //*****************************************************************************************************************************
    GetPythia(topNode, pythia_parton_list, pythia_final_particle_list);

    //*****************************************************************************************************************************
    //Reconstruct Jet
    //*****************************************************************************************************************************
    GetAntiKtParton(pythia_parton_list, pythia_parton_jet_list);
    for(int t = 0; t < pythia_parton_jet_list.size(); t++)
        {
            hPtTrueJet->Fill(pythia_parton_jet_list[t].pT);
        }


    GetAntiKtParticle(pythia_final_particle_list, pythia_final_particle_jet_list);
    for(int r = 0; r < pythia_final_particle_jet_list.size(); r++)
        {
            hPtRecoJet->Fill(pythia_final_particle_jet_list[r].pT);
        }

    //*****************************************************************************************************************************
    //Matched Jet
    //*****************************************************************************************************************************
    GetMatchedJet(pythia_parton_jet_list, pythia_final_particle_jet_list, pythia_matched_jet_list,
                  hDistance, hPtMatchedJet, hPtForEfficiency, hResponseMatrix,
                  hForCorrectionFactor);

    //*****************************************************************************************************************************
    //Fill Trees
    //*****************************************************************************************************************************
    if(false){
    pythiaParton->Fill();
    pythiaFinalParticle->Fill();

    pythiaPartonJet->Fill();
    pythiaFinalParticleJet->Fill();

    pythiaMatchedJet->Fill();
    }

    // any other return code might lead to aborting the event or analysis
    return 0;
}

int JetNloCorrection::End(PHCompositeNode *topNode)
{
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Total events processed = " << nRunEvents << endl;

    if(nRunEvents > 0)
        {
            cout << "This is the end. If everthing was ok, you will see this message: " << endl;
            cout << "123ALLDONE321" << endl;
        }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    hEvents->SetBinContent(1, nRunEvents);

    outfile->Write();
    outfile->Close();

    return 0;
}



void JetNloCorrection::GetPythia(PHCompositeNode *topNode,
                                 std::vector<particles>& parton_list,
                                 std::vector<particles>& final_list)
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

            float energy = part->GetEnergy();
            TLorentzVector *fVector = new TLorentzVector(part->GetPx(), part->GetPy(), part->GetPz(), energy);
            float eta   = fVector->Eta();
            float theta = fVector->Theta();
            float mom   = fVector->P();
            float pT    = mom * sin(theta);
            float eT    = energy * sin(theta);
            float phi   = jetAnalyzer->phiReduce(fVector->Phi());

            //KF- falvor
            int kf = part->GetKF();

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

            //Store the partons before final state particle cuts
            bool passParton = ((phpythia->getLineNumber(part) > 8) && ((part->GetParent() == 7) || (part->GetParent() == 8)));
            //bool passParton = ((phpythia->getLineNumber(part) > 8) && (part->GetParent() <= 8));

            if(passParton)
                {
                    parton_list.push_back(temp);
                }


            //Ignore all beam remnants
            if ( sqrt(pow(part->GetPx(), 2) + pow(part->GetPy(), 2)) == 0.0)
                {
                    continue;
                }
            //Final state particles
            bool passKS = part->GetKS() == 1;

            //No neutrino
            if ((kf == PY_NU_E) || (kf == -PY_NU_E) ||
		(kf == PY_NU_MU) || (kf == -PY_NU_MU) ||
		(kf == PY_NU_TAU) || (kf == -PY_NU_TAU))
                {
                    continue;
                }

            //No muons
            if ((kf == PY_MU) || (kf == -PY_MU))
                {
                    continue;
                }

            //Central acceptance (true = 2pi, false = PHENIX)
            TAcceptParticle *tacceptparticle = new TAcceptParticle();
            tacceptparticle->SetFullTwoPi(false);

            bool inPhi, inEta, inBigEta, deadArea;
            int ok = tacceptparticle->acceptParticle(fVector, charge, zvertex, inPhi, inEta, inBigEta, deadArea);
            delete fVector;

            bool good = ((ok == 1) || (ok == 2));

            //Minimum pT/energy cut on particles
            bool passPtEnergy = false;
            if(charge == 0)
                {
                    passPtEnergy = (energy > CLUSTER_MIN_ENERGY_CUT);
                }
            else
                {
                    passPtEnergy = (pT > TRACK_MIN_PT_CUT);
                }

            //Final state particles in PHENIX acceptance
            if (passKS && good && passPtEnergy)
                {
                    final_list.push_back(temp);
                }
        }
}

void JetNloCorrection::GetAntiKtParton(std::vector<particles> particle_list,
                                       std::vector<jets>& antikt_jet_list)
{
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.2, fastjet::E_scheme, fastjet::Best);

    //Generate input for jet reconstuction
    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    for (unsigned int h = 0; h < particle_list.size(); h++)
        {
            fastjet::PseudoJet pseudoCharged(particle_list[h].px,
                                             particle_list[h].py,
                                             particle_list[h].pz,
                                             particle_list[h].energy);

            jetParticles_all.push_back(pseudoCharged);
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float jetPt  = aFastJet.perp();
            float jetEta = aFastJet.pseudorapidity();
            float jetPhi = jetAnalyzer->phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();

            int arm = 1;
            if (jetPhi > 1.57)
                {
                    arm = 0;
                }

            //if ((fabs(jetEta) < 0.35) && ((jetPhi>-0.58 && jetPhi<0.98) || (jetPhi>2.15 && jetPhi<3.73)))
            if (fabs(jetEta) < 0.35)
                {
                    jets temp;
                    temp.arm           = arm;
                    temp.centralityBin = -99.9;
                    temp.pT            = jetPt;
                    temp.eta           = jetEta;
                    temp.phi           = jetPhi;
                    temp.nc            = (float)nconst;
                    temp.cf            = -99.9;
                    temp.nf            = -99.9;
                    temp.disc          = -99.39;

                    antikt_jet_list.push_back(temp);
                }
        }
    delete antikt;
}

void JetNloCorrection::GetAntiKtParticle(std::vector<particles> particle_list,
					 std::vector<jets>& antikt_jet_list)
{

    std::vector<particles> constituent_list;
    constituent_list.clear();
    jetAnalyzer->GetAntiKtCommon(particle_list, 0.2, 3.0, 10.0, 0.0, 1.0, antikt_jet_list, constituent_list, false);

}

void JetNloCorrection::GetMatchedJet(std::vector<jets> true_jet_list,
                                     std::vector<jets> good_jet_list,
                                     std::vector<jets>& matched_jet_list,
                                     TH1F *hDr,
                                     TH1F *hMatched,
                                     TH1F *hEfficiency,
                                     TH2F *hMatrix,
                                     TH1F *hForCF[17])
{
    std::vector<jetPair> jet_pair_list;
    jet_pair_list.clear();
    unsigned int idTrue = 0;
    for(unsigned int t = 0; t < true_jet_list.size(); t++)
        {
            idTrue++;
            int tArm         = true_jet_list[t].arm;
            float tPt        = true_jet_list[t].pT;
            float tEta       = true_jet_list[t].eta;
            float tPhi       = true_jet_list[t].phi;
            float tNc        = true_jet_list[t].nc;
            float tCf        = true_jet_list[t].cf;
            float tNf        = true_jet_list[t].nf;
            float tDisc      = true_jet_list[t].disc;

            unsigned int idGood = 0;
            for(unsigned int g = 0; g < good_jet_list.size(); g++)
                {
                    idGood++;
                    float gArm       = good_jet_list[g].arm;
                    float gPt        = good_jet_list[g].pT;
                    float gEta       = good_jet_list[g].eta;
                    float gPhi       = good_jet_list[g].phi;
                    float gNc        = good_jet_list[g].nc;
                    float gCf        = good_jet_list[g].cf;
                    float gNf        = good_jet_list[g].nf;
                    float gDisc      = good_jet_list[g].disc;

                    if (tArm == gArm)
                        {
                            float deltaR = jetAnalyzer->dR(tEta, gEta, tPhi, gPhi);

                            jetPair tempPair;
                            tempPair.id.first    = idTrue;
                            tempPair.id.second   = idGood;
                            tempPair.pT.first    = tPt;
                            tempPair.pT.second   = gPt;
                            tempPair.eta.first   = tEta;
                            tempPair.eta.second  = gEta;
                            tempPair.phi.first   = tPhi;
                            tempPair.phi.second  = gPhi;
                            tempPair.nc.first    = tNc;
                            tempPair.nc.second   = gNc;
                            tempPair.cf.first    = tCf;
                            tempPair.cf.second   = gCf;
                            tempPair.nf.first    = tNf;
                            tempPair.nf.second   = gNf;
                            tempPair.disc.first  = tDisc;
                            tempPair.disc.second = gDisc;
                            tempPair.arm         = tArm;
                            tempPair.deltaR      = deltaR;

                            jet_pair_list.push_back(tempPair);
                        }
                }
        }

    //Sort pair by ascending order of deltaR
    std::sort(jet_pair_list.begin(), jet_pair_list.end(), sortPair());

    //Require 1 to 1 matching- and save as unique pair
    std::vector<jetPair> unique_pair_list;
    while (jet_pair_list.size())
        {
            unique_pair_list.push_back(jet_pair_list.front());
            jet_pair_list.erase(std::remove_if(jet_pair_list.begin(), jet_pair_list.end(), removePairId(jet_pair_list.front())),
                                jet_pair_list.end());
        }

    for(unsigned int u = 0; u < unique_pair_list.size(); u++)
        {
            float truePt = unique_pair_list[u].pT.first;
            float recoPt = unique_pair_list[u].pT.second;
            float deltaR = unique_pair_list[u].deltaR;

            hDr->Fill(deltaR);

            hMatched->Fill(recoPt);
            hEfficiency->Fill(truePt);
            hMatrix->Fill(truePt, recoPt);

            jets temp;
            temp.arm           = unique_pair_list[u].arm;
            temp.centralityBin = -1;
            temp.pT            = unique_pair_list[u].pT.second;
            temp.eta           = unique_pair_list[u].eta.second;
            temp.phi           = unique_pair_list[u].phi.second;
            temp.nc            = unique_pair_list[u].nc.second;
            temp.cf            = unique_pair_list[u].cf.second;
            temp.nf            = unique_pair_list[u].nf.second;
            temp.disc          = unique_pair_list[u].disc.second;

            matched_jet_list.push_back(temp);

            float factor = recoPt / truePt;

            if(truePt > PTBINS_FINAL[0]  && truePt <= PTBINS_FINAL[1])
                {
                    hForCF[0]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[1]  && truePt <= PTBINS_FINAL[2])
                {
                    hForCF[1]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[2]  && truePt <= PTBINS_FINAL[3])
                {
                    hForCF[2]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[3]  && truePt <= PTBINS_FINAL[4])
                {
                    hForCF[3]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[4]  && truePt <= PTBINS_FINAL[5])
                {
                    hForCF[4]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[5]  && truePt <= PTBINS_FINAL[6])
                {
                    hForCF[5]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[6]  && truePt <= PTBINS_FINAL[7])
                {
                    hForCF[6]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[7]  && truePt <= PTBINS_FINAL[8])
                {
                    hForCF[7]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[8]  && truePt <= PTBINS_FINAL[9])
                {
                    hForCF[8]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[9]  && truePt <= PTBINS_FINAL[10])
                {
                    hForCF[9]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[10]  && truePt <= PTBINS_FINAL[11])
                {
                    hForCF[10]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[11]  && truePt <= PTBINS_FINAL[12])
                {
                    hForCF[11]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[12]  && truePt <= PTBINS_FINAL[13])
                {
                    hForCF[12]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[13]  && truePt <= PTBINS_FINAL[14])
                {
                    hForCF[13]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[14]  && truePt <= PTBINS_FINAL[15])
                {
                    hForCF[14]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[15]  && truePt <= PTBINS_FINAL[16])
                {
                    hForCF[15]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[16]  && truePt <= PTBINS_FINAL[17])
                {
                    hForCF[16]->Fill(factor);
                }
            if(truePt > PTBINS_FINAL[17]  && truePt <= PTBINS_FINAL[18])
                {
                    hForCF[17]->Fill(factor);
                }

	    if(false){
            if(recoPt > PTBINS_FINAL[0]  && recoPt <= PTBINS_FINAL[1])
                {
                    hForCF[0]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[1]  && recoPt <= PTBINS_FINAL[2])
                {
                    hForCF[1]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[2]  && recoPt <= PTBINS_FINAL[3])
                {
                    hForCF[2]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[3]  && recoPt <= PTBINS_FINAL[4])
                {
                    hForCF[3]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[4]  && recoPt <= PTBINS_FINAL[5])
                {
                    hForCF[4]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[5]  && recoPt <= PTBINS_FINAL[6])
                {
                    hForCF[5]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[6]  && recoPt <= PTBINS_FINAL[7])
                {
                    hForCF[6]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[7]  && recoPt <= PTBINS_FINAL[8])
                {
                    hForCF[7]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[8]  && recoPt <= PTBINS_FINAL[9])
                {
                    hForCF[8]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[9]  && recoPt <= PTBINS_FINAL[10])
                {
                    hForCF[9]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[10]  && recoPt <= PTBINS_FINAL[11])
                {
                    hForCF[10]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[11]  && recoPt <= PTBINS_FINAL[12])
                {
                    hForCF[11]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[12]  && recoPt <= PTBINS_FINAL[13])
                {
                    hForCF[12]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[13]  && recoPt <= PTBINS_FINAL[14])
                {
                    hForCF[13]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[14]  && recoPt <= PTBINS_FINAL[15])
                {
                    hForCF[14]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[15]  && recoPt <= PTBINS_FINAL[16])
                {
                    hForCF[15]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[16]  && recoPt <= PTBINS_FINAL[17])
                {
                    hForCF[16]->Fill(factor);
                }
            if(recoPt > PTBINS_FINAL[17]  && recoPt <= PTBINS_FINAL[18])
                {
                    hForCF[17]->Fill(factor);
                }
	}
        }
}

















