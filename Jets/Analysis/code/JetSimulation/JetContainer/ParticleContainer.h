#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include <phool.h>
#include <PHObject.h>
#include <vector>

typedef struct
{
    int arm;
    int charge;
    float energy;
    float mom;
    float pT;
    float eT;
    float px;
    float py;
    float pz;
    float eta;
    float phi;
    float phiDC;
} iParticles;

class ParticleContainer : public PHObject
{
 public:
    ParticleContainer();
    ~ParticleContainer();
    void Reset();

    void AddEntry(const iParticles &particles)
    {
        particle_list.push_back(particles);
    }
    unsigned int getSize() const
    {
        return particle_list.size();
    }
    int getArm(unsigned int i) const
    {
        return particle_list[i].arm;
    }
    int getCharge(unsigned int i) const
    {
        return particle_list[i].charge;
    }
    float getEnergy(unsigned int i) const
    {
        return particle_list[i].energy;
    }
    float getMom(unsigned int i) const 
    {
        return particle_list[i].mom;
    }
    float getPt(unsigned int i) const 
    {
        return particle_list[i].pT;
    }
    float getEt(unsigned int i) const
    {
        return particle_list[i].eT;
    }
    float getPx(unsigned int i) const
    {
        return particle_list[i].px;
    }
    float getPy(unsigned int i) const
    {
        return particle_list[i].py;
    }
    float getPz(unsigned int i) const
    {
        return particle_list[i].pz;
    }
    float getEta(unsigned int i) const
    {
        return particle_list[i].eta;
    }
    float getPhi(unsigned int i) const
    {
        return particle_list[i].phi;
    }
    float getPhiDC(unsigned int i) const
    {
        return particle_list[i].phiDC;
    }

 private:
    mutable std::vector<iParticles> particle_list;

    ClassDef(ParticleContainer, 1);
};


#endif
