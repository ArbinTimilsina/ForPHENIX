#include <JetContainer.h>

//ClassImp(JetContainer)

JetContainer::JetContainer()
{
    Reset();
}

JetContainer::~JetContainer()
{
    Reset();
}

void JetContainer::Reset()
{
    jet_list.clear();
}


