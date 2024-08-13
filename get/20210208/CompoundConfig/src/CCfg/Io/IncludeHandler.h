// =============================================================================
// IncludeHandler.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Io_IncludeHandler_h_INCLUDED
#define CCfg_Io_IncludeHandler_h_INCLUDED

#include <string>
#include <memory>

namespace CCfg
{
// forward declarations
class Object;

namespace Io
{
// forward declarations
class Url;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * A callback interface to take custom actions when an include is needed while
 * parsing or formatting.
 * This callback is used both when a parser encounters an include directive in
 * a configuration, and when a formatter encounters a non empty include URL in
 * a configuration tree.
 *
 * @author fchateau
 */
class IncludeHandler
{
public:
	// construction/destruction
	IncludeHandler();
	virtual ~IncludeHandler();

	// callback methods
	/**
	 * Callback method that must be implemented to react include directive
	 * in files.
	 * The object doing the include is created, initialized with the URL and
	 * passed to this callback. Its is up to the implementation to load the
	 * included configuration and attach it to the object.
	 * @param includerObject the object with an include directive. Its name,
	 * index, and include URL are defined.
	 */
	virtual void parseInclude(Object& includerObject) = 0;

	/**
	 * Callback method that must be implemented to react include URL when
	 * formatting a tree structure
	 * The object having a non empty include property is passed to this
	 * callback. The object is the root of the subtree that will have to be
	 * saved in the appropriate place.
	 * @param includerObject the object with an include property.
	 * @return Returns true if the include was processed, false otherwise.
	 * At typical implementation of a formatter would stop recursing in the
	 * includer object if the function returns true.
	 * On the other hand, a false value means that the object content
	 * should be written as if there were no include.
	 */
	virtual bool formatInclude(const Object& includerObject) = 0;
};

/**
 * This is an implementation of IncludeHandler that does purposely nothing.
 * It ignores all include directives.
 *
 * @author fchateau
 */
class NoOpIncludeHandler : public IncludeHandler
{
public:
	~NoOpIncludeHandler();

	// callback methods
	void parseInclude(Object&);
	bool formatInclude(const Object&);

	static NoOpIncludeHandler& getInstance();

private:
	NoOpIncludeHandler();
	static std::auto_ptr<NoOpIncludeHandler> instance;
};

}

}

#endif // CCfg_Io_IncludeHandler_h_INCLUDED
