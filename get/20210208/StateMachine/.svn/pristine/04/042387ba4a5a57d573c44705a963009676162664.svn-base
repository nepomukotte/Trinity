// =============================================================================
// Engine.cpp
// -----------------------------------------------------------------------------
// class dhsm::Engine
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "QueuedEngine.h"
#include "State.h"
#include "utl/Logging.h"

namespace dhsm
{
// =============================================================================
// dhsm::Engine class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
QueuedEngine::QueuedEngine()
	: eventQueue(), stepCount(0)
{
}

/**
 * Destructor.
 */
QueuedEngine::~QueuedEngine()
{
}

/*
 * (see inherited doc)
 */
QueuedEngine::ProcessingResult QueuedEngine::start()
{
	stepCount = 0;
	eventQueue.clear(); // Emptying event queue in case there are remains from a previous run.
	return Engine::start();
}

/*
 * (see inherited doc)
 */
QueuedEngine::ProcessingResult QueuedEngine::step(const Event& evt)
{
	broadcast(evt.clone());
	stepAll();
	return Processed;
}

/**
 * Broadcasts the specified Event into the engine's queue.
 * It will be handled when a future call to step() consumes enough event to
 * reach this new one.
 */
void QueuedEngine::broadcast(EventPtr event)
{
	if(isInitiated(*stateMachine_))
	{
		LOG_DEBUG() << "BCST: " << *event;
		eventQueue.push_front(event.release());
	}
	else
	{
		LOG_ERROR() << "BCST: " << "failed because state machine is not initiated.";
	}
}

/**
 * Broadcast the event with the specified name into the engine's queue.
 * @see void Engine::broadcast(EventPtr&)
 */
void QueuedEngine::broadcast(const char* evtName)
{
	EventPtr event(new Event(evtName));
	broadcast(event);
}

/**
 * Consumes up to the specified number of events in the event queue and triggers all
 * associated Transitions.
 * @param nSteps specifies the number of steps to be executed. A zero or negative
 * number means an unlimited number of steps (until event queue is empty).
 * @return Actual number of Events processed: necessarily smaller than or equal
 * to one specified.
 */
int QueuedEngine::step(int nSteps)
{
	int processedSteps = 0;
	EventPtr event;
	while((event.reset(nextEvent().release()), event.get() != 0) and (nSteps == 0 or processedSteps < nSteps))
	{
		Engine::step(*event);
		++processedSteps;
	}
	return processedSteps;
}

/**
 * @fn int Engine::stepAll()
 * Consumes all Events until Event queue is empty and triggers all associated
 * Transitions.
 * @return Actual number of Events processed.
 */

/**
 * @fn Engine::getStepCount() const
 * Returns the number of steps processed since the engine was enabled.
 */

#ifdef UTL_LOGGING_ENABLED
/**
 * Printouts info on state machine.
 */
void QueuedEngine::dump() const
{
}
#endif

/**
 * Pops the Event queue.
 * @return Returns the popped Event or an empty Event if the queue was empty.
 */
EventPtr QueuedEngine::nextEvent()
{
	EventPtr evt;
	if(not eventQueue.empty())
	{
		evt.reset(eventQueue.pop_front().release());
		stepCount++;
		LOG_DEBUG() << "---> Popped: " << evt->label();
	}
	return evt;
}

} // namespace dhsm
