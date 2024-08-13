// =============================================================================
// Event.h
// -----------------------------------------------------------------------------
// class dhsm::Event
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Event_h_INCLUDED
#define dhsm_Event_h_INCLUDED

#include <ostream>
#include <string>
#include <memory>

namespace dhsm
{
class Event;
typedef std::auto_ptr<Event> EventPtr;

/**
 * This class represent the events: they stand for some condition that must be
 * fulfilled to cross a transition.
 * An Event is characterized by its label. Different events should not have the
 * same label.
 * Events are issued to the state machine's event queue by calls to the
 * Engine::broadcast() method, either within a transition or directly. Events
 * are consumed by the state machine's engine whenever the Engine::step() method
 * is called.
 */
class Event
{
public:
	Event();
	Event(const char* label);
	virtual ~Event();

	virtual EventPtr clone() const;

	bool isValid() const { return not label().empty(); }

	bool operator< (const Event& event) const { return label() <  event.label(); }
	bool operator==(const Event& event) const { return label() == event.label(); }

	const std::string& label() const { return label_; }
	void label(const std::string& newLabel) { label_ = newLabel; }

private:
	std::string label_;
};

inline std::ostream& operator<<(std::ostream& stream, const Event& event)
{
	stream << '#' << event.label() << '#';
	return stream;
}

} // namespace dhsm

#endif // dhsm_Event_h_INCLUDED
