#ifndef JETCONTAINER_H
#define JETCONTAINER_H

#include <phool.h>
#include <PHObject.h>
#include <vector>

typedef struct
{
    int arm;
    float pT;
    float eta;
    float phi;
    float nc;
    float cf;
    float disc;
} iJets;

class JetContainer : public PHObject
{
 public:
    JetContainer();
    ~JetContainer();
    void Reset();

    void AddEntry(const iJets jets)
    {
        jet_list.push_back(jets);
    }
    unsigned int getSize() const
    {
        return jet_list.size();
    }
    int getArm(unsigned int i) const
    {
        return jet_list[i].arm;
    }
    float getPt(unsigned int i) const
    {
        return jet_list[i].pT;
    }
    float getEta(unsigned int i) const
    {
        return jet_list[i].eta;
    }
    float getPhi(unsigned int i) const
    {
        return jet_list[i].phi;
    }
    float getNc(unsigned int i) const
    {
        return jet_list[i].nc;
    }
    float getCf(unsigned int i) const
    {
        return jet_list[i].cf;
    }
    float getDisc(unsigned int i) const
    {
        return jet_list[i].disc;
    }

 private:
    mutable std::vector<iJets> jet_list;

    ClassDef(JetContainer, 1);
};


#endif
