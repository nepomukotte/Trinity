// =============================================================================
// Format.h
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

#ifndef CCfg_Io_Format_h_INCLUDED
#define CCfg_Io_Format_h_INCLUDED

#include "../Exception.h"
#include <memory>
#include <string>

namespace CCfg
{
// forward declarations
class Object;

namespace Io
{
// foward declarations
class Storage;
class IncludeHandler;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * An interface that each configuration serialization format must realize.
 *
 * @author fchateau
 */
class Format
{
public:
	// inner classes
	class InvalidFormat;
	class SyntaxError;
	class SemanticError;

	// construction/destruction
	Format();
	virtual ~Format();

	std::auto_ptr<CCfg::Object> parse(Storage& storage);
	void format(const Object *config, Storage& storage);

	virtual std::auto_ptr<CCfg::Object> parse(Storage& storage, IncludeHandler& handler) = 0;
	virtual void format(const Object *config, Storage& storage, IncludeHandler& handler) = 0;
	virtual std::string getConfigName(Storage& storage) = 0;
};

/**
 * Reports that a document cannot be parsed because it does not respect
 * format specifications.
 */
class Format::InvalidFormat : public Exception
{
public:
	// construction/destruction
	InvalidFormat();
	virtual ~InvalidFormat() throw();
};

/**
 * Reports that a document cannot be parsed because there is a syntax error.
 */
class Format::SyntaxError : public Format::InvalidFormat
{
public:
	// construction/destruction
	SyntaxError();
	virtual ~SyntaxError() throw();
};

/**
 * Reports that a document cannot be parsed because the semantics constraints
 * on the format are not respected.
 */
class Format::SemanticError : public Format::InvalidFormat
{
public:
	// construction/destruction
	SemanticError();
	virtual ~SemanticError() throw();
};

}
}

#endif // CCfg_Io_Format_h_INCLUDED
