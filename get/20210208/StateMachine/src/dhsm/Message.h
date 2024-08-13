// =============================================================================
// Message.h
// -----------------------------------------------------------------------------
// class dhsm::Message
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Message_h_INCLUDED
#define dhsm_Message_h_INCLUDED

#include "Event.h"

namespace dhsm
{
/**
 * This class extends the events to allow passing one parameter alongside with
 * the event.
 * The callback method for transition actions and enter actions that may be
 * triggered by the event must have a parameter of the same type in their
 * signature to be able to receive the value.
 * It is also perfectly valid for some callbacks to have a void signature to
 * ignore the parameter.
 * However, if the callback excepts an argument of another type, a runtime_error
 * will be raised.
 */
template<typename T>
class Message : public Event
{
	T parameter_;
public:
	Message() {}
	Message(const char* eventName, const T& parameter) : Event(eventName), parameter_(parameter) {}
	~Message() {}

	EventPtr clone() const { return EventPtr(new Message<T>(*this)); }

	T& parameter() { return parameter_; }
	const T& parameter() const { return parameter_; }
};

/**
 * Instanciates a new Message on the heap, initializes it and return it.
 */
template<typename T>
EventPtr createMessage(const char* eventName, const T& parameter)
{
	return EventPtr(new Message<T>(eventName, parameter));
}

}

#endif // dhsm_Message_h_INCLUDED
