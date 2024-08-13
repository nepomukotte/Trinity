// =============================================================================
// QueuedEngine.h
// -----------------------------------------------------------------------------
// class dhsm::QueuedEngine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_QueuedEngine_h_INCLUDED
#define dhsm_QueuedEngine_h_INCLUDED

#include "Engine.h"
#include <boost/ptr_container/ptr_deque.hpp>
#include <string>

namespace dhsm
{
/**
 * This class extends the Engine by adding an event queue to support
 * asynchronous event processing.
 *
 * Any engine instance has its own queue of events. By calling its step()
 * method, you consume an event and trigger the corresponding Transition.
 * The engine then executes the list of actions associated to the transition
 * then enters the target state specified in the transition.
 * By calling its broadcast() method you feed the event queue with a new Event.
 *
 * @author sanvar, fchateau.
 */
class QueuedEngine : public Engine
{
public:
	// Typedefs
	typedef boost::ptr_deque<Event> EventQueue;

public:
	QueuedEngine();
	virtual ~QueuedEngine();
	//@{
	void broadcast(EventPtr event);
	void broadcast(const char* evtName);
	virtual ProcessingResult start();
	virtual ProcessingResult step(const Event& event);
	//@}
	//@{
	int step(int nSteps = 1);
	int stepAll() { return step(0); }
	int getStepCount() const { return stepCount; }
	//@}

	#ifdef UTL_LOGGING_ENABLED
	void dump() const;
	#endif

private:
	EventPtr nextEvent();

private:
	EventQueue eventQueue; ///< The queue of events waiting to be processed.
	int stepCount; ///< The number of processed steps since the engine was initiated.
};

} // namespace dhsm

#endif // dhsm_QueuedEngine_h_INCLUDED
