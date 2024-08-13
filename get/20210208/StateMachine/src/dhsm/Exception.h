// =============================================================================
// Exception.h
// -----------------------------------------------------------------------------
// class dhsm::Exception
// class dhsm::StaticException
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Exception_h_INCLUDED
#define dhsm_Exception_h_INCLUDED

#include <utl/Exception.h>

namespace dhsm
{
/**
 * This is the base class of all exceptions raised by the state machine framework.
 */
class Exception : public utl::Exception
{
public:
	Exception();
	virtual ~Exception() throw();

	virtual void writeMessage(std::ostream& stream) const = 0;
};

/**
 * Exceptions heriting from this one are traversing the Engine::process handlers
 * without being caught.
 */
class EngineTraversalException : public Exception
{
public:
	EngineTraversalException();
	virtual ~EngineTraversalException() throw();
};

/**
 * This exception signals an error that cannot be handled internally by the
 * state machine error mecanism.
 * All exceptions that should pass through the framework handlers and fall
 * back to the caller should subclass it.
 */
class CriticalError : public EngineTraversalException
{
public:
	CriticalError();
	virtual ~CriticalError() throw();

	virtual void writeMessage(std::ostream& stream) const = 0;
};

/**
 * This exception is thrown internally to signal an error during a transition
 * that have been recovered.
 * Its exists mainly to make the exception cross the Engine::process handler stack
 * all the way down to Engine::step method, once the initial std::exception has
 * been catched.
 */
class RecoveredError : public EngineTraversalException
{
public:
	RecoveredError(const std::string& msg);
	virtual ~RecoveredError() throw();

	virtual void writeMessage(std::ostream& stream) const;

private:
	std::string message;
};

/**
 * This exception signals an error inherent to the declaration of a state machine.
 */
class StaticError : public CriticalError
{
public:
	enum Kind { MessageTypeMismatch,
		ClosureClassMismatch,
		MissingInitialTransition,
		UnsupportedEvent,
		NoSuchState,
		NoSuchStateMachine,
		NotCompositeState
	};

	StaticError(Kind errorKind, const std::string& detail);
	~StaticError() throw();

	void writeMessage(std::ostream& stream) const;

protected:
	Kind errorKind;
	std::string detail;
};

class TransitionError : public Exception
{
public:
	explicit TransitionError(const std::string& msg);
	virtual ~TransitionError() throw();

	virtual void writeMessage(std::ostream& stream) const;

private:
	std::string message;
};

}

#endif // dhsm_Exception_h_INCLUDED
