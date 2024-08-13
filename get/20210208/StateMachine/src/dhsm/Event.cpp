// =============================================================================
// Event.cpp
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

#include "Event.h"

namespace dhsm
{
// =============================================================================
// dhsm::Event class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
Event::Event() : label_()
{
}

/**
 * Full constructor with the event name.
 */
Event::Event(const char* label) : label_(label)
{
}

Event::~Event()
{
}

EventPtr Event::clone() const
{
	return EventPtr(new Event(*this));
}

/**
 * @fn bool Event::isValid() const
 * Returns whether this event is valid (ie initialized with a non-empty label).
 */

/**
 * @fn bool Event::operator< (const Event& event) const
 * Provide an ordering according to the lexicographical order of their labels.
 */

/**
 * @fn bool Event::operator==(const Event& event) const
 * Check the equality of two events according to the equality of their labels.
 */

} // namespace dhsm
