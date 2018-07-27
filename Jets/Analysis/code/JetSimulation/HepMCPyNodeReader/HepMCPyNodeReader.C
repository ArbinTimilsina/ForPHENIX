//Thanks to Mike McCumber for this
//Some modifications made by Arbin
#include "HepMCPyNodeReader.h"
//#include "PHG4InEvent.h"
//#include "PHG4Particlev1.h"

#include <PHHepMCGenEvent.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,8)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

#include <PHHijingHeader.h>
#include <PHHijingHeaderV3.h>
#include <PHPythiaContainer.h>
#include <PHPythiaContainerV2.h>

#include <PHIODataNode.h>
#include <PHObject.h>
#include <PHCompositeNode.h>
#include <PHNodeIterator.h>
#include <PHNodeReset.h>
#include <PHTimeStamp.h>
#include <Fun4AllReturnCodes.h>

#include <Fun4AllReturnCodes.h>
#include <getClass.h>
#include <recoConsts.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <gsl/gsl_const.h>

#include <list>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

HepMCPyNodeReader::HepMCPyNodeReader(const std::string &name):
    SubsysReco(name)
{
    hijingHeader = NULL;
    pythiaContainer = NULL;
}

int HepMCPyNodeReader::InitRun(PHCompositeNode *topNode)
{

    cout << "InitRun called in HepMCPyNodeReader" << endl;

    //Initialize the PHPythia container
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
        {
            cout << PHWHERE << "DST Node missing doing nothing" << endl;
            return ABORTRUN;
        }

    //Hijing header information
    hijingHeader = new PHHijingHeaderV3();
    PHObjectNode_t *PHHijingHeaderNode = new PHObjectNode_t(hijingHeader, "PHHijingHeader", "PHObject");
    dstNode->addNode(PHHijingHeaderNode);

    //Hijing/Pythia particle information
    pythiaContainer = new PHPythiaContainerV2();
    PHObjectNode_t *PHPythiaNode = new PHObjectNode_t(pythiaContainer, "PHHijing", "PHObject");
    dstNode->addNode(PHPythiaNode);

    return 0;
}

int HepMCPyNodeReader::process_event(PHCompositeNode *topNode)
{
    //Find the HepMC data
    PHHepMCGenEvent *genEvt = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
    HepMC::GenEvent *hepMC =  genEvt->getEvent();
    if (!hepMC)
        {
            cout << PHWHERE << " no hepMC pointer under HEPMC Node found" << endl;
            return ABORTRUN;
        }

    //In HepMC default value is set by configure switch, so convert to the units we want
    //TMCParticle has Energy [GeV], Momenta [GeV/c], and Vertex [mm]
    const double momFactor = HepMC::Units::conversion_factor( hepMC->momentum_unit(), HepMC::Units::GEV );
    const double lengthFactor = HepMC::Units::conversion_factor( hepMC->length_unit(), HepMC::Units::MM );

    //In HepMC Each vertex maintains a listing of its incoming and outgoing particles,
    //while each particle points back to its production vertex and decay vertex.

    //Use the first vertex encountered as the primary
    double xPos = 0.0;
    double yPos = 0.0;
    double zPos = 0.0;
    bool first = true;

    unsigned int nparticles = 0;
    //GenEvent::vertex_iterator: It walks through all vertices in the event exactly once
    for ( HepMC::GenEvent::vertex_iterator v = hepMC->vertices_begin(); v != hepMC->vertices_end(); ++v )
        {
            if (first)
                {
                    xPos = (*v)->position().x() * lengthFactor;
                    yPos = (*v)->position().y() * lengthFactor;
                    zPos = (*v)->position().z() * lengthFactor;

                    first = false;
                }
            //GenVertex::particle_iterator- children: Walks over all particles outgoing from the root/vertex
            for (HepMC::GenVertex::particle_iterator p = (*v)->particles_begin(HepMC::children); p != (*v)->particles_end(HepMC::children); ++p)
                {
                    const HepMC::FourVector& momVector = (*p)->momentum();

                    //We want to pass only the final state particles through PISA
                    //1. Don't keep the particle if there is decay vertex
                    //2. Don't keep the particle that has decayed (status 2) + or is Fragmentation Jet (status 103)
                    if (((*p)->end_vertex()) || ((*p)->status() != 1))
                        {
                            continue;
                        }

		    float particleRho = sqrt(pow(momVector.px() * momFactor, 2) +
                                             pow(momVector.py() * momFactor, 2) +
                                             pow(momVector.pz() * momFactor, 2));
                    float particleEta = 0.5 * log((particleRho + (momVector.pz() * momFactor)) / (particleRho - (momVector.pz() * momFactor)));

		    //Save time while passing through PISA- don't care about Muon arm or BBC/ZDC
                    if(fabs(particleEta) > 1.0)
                        {
                            continue;
                        }

                    //Set the particle information
                    TMCParticle particle;
                    particle.SetEnergy(momVector.e()*momFactor);
                    particle.SetKF((*p)->pdg_id());
                    particle.SetKS((*p)->status());
                    particle.SetMass(momVector.m()*momFactor);
                    particle.SetPx(momVector.px()*momFactor);
                    particle.SetPy(momVector.py()*momFactor);
                    particle.SetPz(momVector.pz()*momFactor);
                    particle.SetVx((*v)->position().x()*lengthFactor);
                    particle.SetVy((*v)->position().y()*lengthFactor);
                    particle.SetVz((*v)->position().z()*lengthFactor);

                    particle.SetParent(-1);
                    particle.SetFirstChild(-1);
                    particle.SetLastChild(-1);
                    particle.SetTime(0.0);

                    pythiaContainer->addParticle(particle);
                    ++nparticles;
                }
        }

    //Write the header information
    HepMC::HeavyIon* heavyIon = hepMC->heavy_ion();

    hijingHeader->SetEvt(hepMC->event_number());                       //Event number
    hijingHeader->SetNpart(nparticles);                                //Number of particles in event
    hijingHeader->SetATarg(0);                                         //A of target
    hijingHeader->SetZTarg(0);                                         //Z of target
    hijingHeader->SetAProj(0);                                         //Projectile A
    hijingHeader->SetZProj(0);                                         //Projectile Z
    hijingHeader->SetCollisionE(0);                                    //sqrt(s_nn) of collision
    hijingHeader->SetNbinary(heavyIon->Ncoll());                       //Number of binary collisions
    hijingHeader->SetBimpact(heavyIon->impact_parameter());            //Impact Parameter
    hijingHeader->SetNp(heavyIon->Npart_proj());                       //Number of participating nucleons in projectile
    hijingHeader->SetNt(heavyIon->Npart_targ());                       //Number of participating nucleons in target

    hijingHeader->SetPrimaryVertexX(xPos);
    hijingHeader->SetPrimaryVertexY(yPos);
    hijingHeader->SetPrimaryVertexZ(zPos);

    return EVENT_OK;
}

int HepMCPyNodeReader::ResetEvent(PHCompositeNode *topNode)
{
    PHNodeIterator mainIter(topNode);
    PHNodeReset reset;
    if (mainIter.cd(Name()))
        {
            mainIter.forEach(reset);
        }
    return 0;
}












