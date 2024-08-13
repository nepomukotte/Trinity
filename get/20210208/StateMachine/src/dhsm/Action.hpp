// =============================================================================
// Action.hpp
// -----------------------------------------------------------------------------
// class dhsm::AbstractAction
// class dhsm::Action<T>
// class dhsm::ActionList
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_Action_hpp_INCLUDED
#define dhsm_Action_hpp_INCLUDED

#include "StateData.h"
#include "Message.h"
#include "Exception.h"
#include <utl/Logging.h>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <ostream>
#include <string>
#include <stdexcept>
#include <typeinfo>

namespace dhsm
{
/**
 * This is the base class of all actions executed in response to state machine
 * events like transition, stater enter/exit, reset.
 * There are three main kinds of Actions in a state machine:
 * - actions associated to transitions, executed when a transition is triggered
 * - actions associated to states, executed when a state is entered or left,
 * - actions associated to state machines, executed when it is reset.
 * AbstractAction provides a virtual execute() method that all deriving
 * classes are expected to override. This method is meant to encapsulate the
 * operations corresponding to an actual Action.
 */
class AbstractAction
{
public:
	virtual ~AbstractAction();
	virtual void execute(const Event& event, StateData* closureData) const = 0;

	#ifdef UTL_LOGGING_ENABLED
	const std::string& label() const { return label_; }
	void label(const std::string& newLabel) { label_ = newLabel; }
	#endif

protected:
	template<class T, class U> struct Method;

private:
	#ifdef UTL_LOGGING_ENABLED
	std::string label_;
	#endif
};

template<class T, class U>
void callAction(T& instance, void (T::*callback)(U), const Event& event);

/**
 * This template class allows users to wrap any object method having the
 * signature 'void (T::*)(U)' in an action.
 * This allows users to avoid going through the hassle of defining one class for
 * every Action they want to implement in their application.
 */
template<class T, class U>
class Action : public AbstractAction
{
public:
	typedef typename Method<T,U>::Signature Callback; ///< Signature of the action method.

public:
	Action(T& instance, Callback callback);
	void execute(const Event& event, StateData* closureData) const;

protected:
	T& instance; ///< Instance on which the callback method is expected to run.
	Callback callback; ///< Callback method associated to this action.
};

/**
 * This is a closure action, the method wrapped inside this class will be
 * executed in the context of the current StateData.
 * The user should be careful about the class where this method belongs:
 *  - enter and exit closures are executed using the state data associated with
 * the state we enter/exit.
 *  - transition and reset actions are executed using the state data of the
 * enclosing composite state.
 */
template<class T, class U>
class Closure : public AbstractAction
{
public:
	typedef typename Method<T,U>::Signature Callback; ///< Signature of the closure method.

public:
	Closure(Callback callback);
	void execute(const Event& event, StateData* closureData) const;

protected:
	Callback callback; ///< Callback method associated to this action.
};

/**
 * This class applies the "composite pattern" to regroup multiple Actions.
 * It shares the same interface as Action and calling ActionList::execute()
 * causes all the Actions it contains to execute (in order of creation).
 *
 * If any action reports a failure, the next actions are not executed.
 */
class ActionList : public AbstractAction
{
public:
	void clear() { actions.clear(); }
	size_t count() const { return actions.size(); }
	ActionList& add(AbstractAction* action);
	template<class T>
		ActionList& add(T& instance, void (T::*callback)());
	template<class T, class U>
		ActionList& add(T& instance, void (T::*callback)(U));
	template<class T>
		ActionList& add(void (T::*callback)());
	template<class T, class U>
		ActionList& add(void (T::*callback)(U));

	void execute(const Event& event, StateData* closureData) const;

	#ifdef UTL_LOGGING_ENABLED
	void dump(int indent = 0) const;
	#endif

private:
	typedef boost::ptr_vector<AbstractAction> ActionVector;
	ActionVector actions;	///< The list of actions to execute.
};

/**
 * This class defines the signature of the action's callback method according
 * to its parameter.
 */
template<class T, class U>
struct AbstractAction::Method
{
	typedef void (T::*Signature)(U);
};

template<class T>
struct AbstractAction::Method<T, void>
{
	typedef void (T::*Signature)();
};


// =============================================================================
// dhsm::AbstractAction inline methods body
// -----------------------------------------------------------------------------

#ifdef UTL_LOGGING_ENABLED
inline std::ostream& operator<<(std::ostream& stream, const AbstractAction& action)
{
	stream << '{' << action.label() << '}';
	return stream;
}
#endif

// =============================================================================
// Template global functions
// -----------------------------------------------------------------------------

/**
 * Calls the action callback, extracting the action parameter from the event if necessary.
 */
template<class T, class U>
void callAction(T& instance, void (T::*callback)(U), const Event& event)
{
	typedef typename boost::remove_const<typename boost::remove_reference<U>::type>::type ArgType;
	const Message<ArgType>* msg = dynamic_cast<const Message<ArgType>*>(&event);
	if(msg == 0)
	{
		//@TODO: send more detailed informations
		throw StaticError(StaticError::MessageTypeMismatch, typeid(event).name());
	}
	(instance.*callback)(msg->parameter());
}

/*
 * Specialization of the previous method, handles the case where the action
 * does not have a parameter.
 */
template<class T>
void callAction(T& instance, void (T::*callback)(), const Event& /*event*/)
{
	(instance.*callback)();
}

// =============================================================================
// dhsm::Action class template methods body
// -----------------------------------------------------------------------------

/**
 * Constructor from a method.
 * @param callback Callback method associated to this action.
 * @param instance The T instance on which the callback is called.
 */
template<class T, class U>
Action<T, U>::Action(T& instance, Callback callback)
	: AbstractAction(), instance(instance), callback(callback)
{
}

/**
 * Executes the action's callback method.
 */
template<class T, class U>
void Action<T, U>::execute(const Event& event, StateData* /*closureData*/) const
{
	callAction(instance, callback, event);
}

// =============================================================================
// dhsm::Closure class template methods body
// -----------------------------------------------------------------------------

/**
 * Constructor.
 */
template<class T, class U>
Closure<T, U>::Closure(Callback callback)
	: callback(callback)
{
}

/**
 * Executes the closure using the specified state data.
 */
template<class T, class U>
void Closure<T, U>::execute(const Event& event, StateData* closureData) const
{
	T* object = dynamic_cast<T*>(closureData);
	if(object == 0)
	{
		throw StaticError(StaticError::ClosureClassMismatch, typeid(closureData).name());
	}
	else
	{
		callAction(*object, callback, event);
	}
}

// =============================================================================
// dhsm::ActionList class template methods body
// -----------------------------------------------------------------------------

/**
 * Creates an action calling the specified method on the specified object, and
 * adds it to the list.
 * @return Returns *this to allow chained calls
 */
template<class T>
ActionList& ActionList::add(T& instance, void (T::*callback)())
{
	add(new Action<T, void>(instance, callback));
	return *this;
}

/**
 * Creates a parameterized action calling the specified method on the specified
 * object, and adds it to the list.
 * @return Returns *this to allow chained calls
 */
template<class T, class U>
ActionList& ActionList::add(T& instance, void (T::*callback)(U))
{
	add(new Action<T, U>(instance, callback));
	return *this;
}

/**
 * Creates a closure calling the specified method, and adds it to the list.
 * @return Returns *this to allow chained calls
 */
template<class T>
ActionList& ActionList::add(void (T::*callback)())
{
	add(new Closure<T, void>(callback));
	return *this;
}

/**
 * Creates a parameterized closure calling the specified method, and adds it to
 * the list.
 * @return Returns *this to allow chained calls
 */
template<class T, class U>
ActionList& ActionList::add(void (T::*callback)(U))
{
	add(new Closure<T, U>(callback));
	return *this;
}

} // namespace dhsm

#endif // dhsm_Action_hpp_INCLUDED
