#include <EventContainerV1.h>

ClassImp(EventContainerV1)

EventContainerV1::EventContainerV1()
{
    Reset();
}

EventContainerV1::~EventContainerV1()
{
    Reset();
}

void EventContainerV1::Reset()
{
    event_list.clear();
}

