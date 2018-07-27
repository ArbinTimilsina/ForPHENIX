#ifndef HEPMCPYNODEREADER_H__
#define HEPMCPYNODEREADER_H__

#include <SubsysReco.h>
#include <PHHijing.h>
#include <PHHijingHeader.h>
#include <PHPythia.h>
#include <PHPythiaContainer.h>

#include <string>

class PHCompositeNode;

class HepMCPyNodeReader : public SubsysReco
{

 public:
    HepMCPyNodeReader(const std::string &name = "HEPMCPYNODEREADER");
    virtual ~HepMCPyNodeReader() {}

    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int ResetEvent(PHCompositeNode *topNode);

 private:

    PHHijingHeader* hijingHeader;
    PHPythiaContainer* pythiaContainer;
};

#endif




