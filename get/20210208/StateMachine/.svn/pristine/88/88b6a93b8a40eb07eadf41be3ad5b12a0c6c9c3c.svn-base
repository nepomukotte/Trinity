// =============================================================================
// Exception.cpp
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

#include "Exception.h"

namespace dhsm
{
// =============================================================================
// dhsm::Exception class methods body
// -----------------------------------------------------------------------------

Exception::Exception()
{
}

Exception::~Exception() throw()
{
}
// =============================================================================
// dhsm::EngineTraversalException class methods body
// -----------------------------------------------------------------------------

EngineTraversalException::EngineTraversalException()
{
}

EngineTraversalException::~EngineTraversalException() throw()
{
}

// =============================================================================
// dhsm::CriticalError class methods body
// -----------------------------------------------------------------------------

CriticalError::CriticalError()
{
}

CriticalError::~CriticalError() throw()
{
}

// =============================================================================
// dhsm::StaticError class methods body
// -----------------------------------------------------------------------------

StaticError::StaticError(Kind errorKind, const std::string& detail) : errorKind(errorKind), detail(detail)
{
}

StaticError::~StaticError() throw()
{
}

void StaticError::writeMessage(std::ostream& stream) const
{
	stream << "Static error: ";
	switch(errorKind)
	{
	case MessageTypeMismatch:stream << "unexpected message type";
		break;
	case ClosureClassMismatch:stream << "invalid closure class";
		break;
	case MissingInitialTransition:stream << "missing initial transition";
		break;
	case UnsupportedEvent:stream << "unsupported event";
		break;
	case NoSuchState:stream << "no such state";
		break;
	case NoSuchStateMachine:stream << "no such state-machine";
		break;
	case NotCompositeState:stream << "not a composite state";
		break;
	}
	stream << ": " << detail;
}

// =============================================================================
// dhsm::RecoveredError class methods body
// -----------------------------------------------------------------------------

RecoveredError::RecoveredError(const std::string& msg)
	: message(msg)
{
}

RecoveredError::~RecoveredError() throw()
{
}

void RecoveredError::writeMessage(std::ostream& stream) const
{
	stream << message;
}

// =============================================================================
// dhsm::TransitionError class methods body
// -----------------------------------------------------------------------------

TransitionError::TransitionError(const std::string& msg)
	: message(msg)
{
}

TransitionError::~TransitionError() throw()
{
}

void TransitionError::writeMessage(std::ostream& stream) const
{
	stream << message;
}

}
