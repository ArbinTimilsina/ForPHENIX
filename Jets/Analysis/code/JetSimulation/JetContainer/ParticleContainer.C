#include <ParticleContainer.h>

//ClassImp(ParticleContainer)

ParticleContainer::ParticleContainer()
{
    Reset();
}

ParticleContainer::~ParticleContainer()
{
    Reset();
}

void ParticleContainer::Reset()
{
    particle_list.clear();
}


