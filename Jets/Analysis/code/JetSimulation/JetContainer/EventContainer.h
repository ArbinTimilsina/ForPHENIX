#ifndef EVENTCONTAINER_H
#define EVENTCONTAINER_H

#include <phool.h>
#include <PHObject.h>
#include <vector>

typedef struct
{
    float centrality;
    float zVertexPythia;
    float zVertexPisaReco;
    float zVertexData;
}iEvents;

class EventContainer : public PHObject
{
 public:
    EventContainer();
    virtual ~EventContainer();

    virtual void Reset()
    { PHOOL_VIRTUAL_WARNING; }

    virtual void AddEntry(const iEvents events)
    {
	PHOOL_VIRTUAL_WARNING; return; 
    }

    virtual unsigned int getSize() const
    {
	PHOOL_VIRTUAL_WARNING; return 0; 
    }

    virtual float getCentrality(unsigned int i) const
    {
	PHOOL_VIRTUAL_WARNING; return 0.0; 
    }
    virtual float getVertexPythia(unsigned int i) const
    {
        PHOOL_VIRTUAL_WARNING; return 0.0;
    }
    virtual float getVertexPisaReco(unsigned int i) const
    {
        PHOOL_VIRTUAL_WARNING; return 0.0;
    }
    virtual float getVertexData(unsigned int i) const
    {
        PHOOL_VIRTUAL_WARNING; return 0.0;
    }

 private:

    ClassDef(EventContainer, 1);
};


#endif
