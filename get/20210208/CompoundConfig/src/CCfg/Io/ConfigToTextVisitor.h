// =============================================================================
// ConfigToTextVisitor.h
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

#ifndef CCfg_Io_ConfigToTextVisitor_h_INCLUDED
#define CCfg_Io_ConfigToTextVisitor_h_INCLUDED

#include "../Visitor.h"
#include <string>
#include <ostream>
#include <vector>

namespace CCfg
{
// forward declarations
class ClassInfo;

namespace Io
{
// forward declarations
class IncludeHandler;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for classes that output a configuration tree to a
 * text stream.
 *
 * @author fchateau
 */
class ConfigToTextVisitor : public Visitor
{
protected:
	typedef std::vector<const ClassInfo*> WrittenClassList;

	IncludeHandler& includeHandler;
	std::ostream& stream;
	std::string indentPattern;

	WrittenClassList writtenClasses;

	// transient attributes
	unsigned indentLevel;
	bool nextIsRootObject;

public:
	// construction/destruction
	ConfigToTextVisitor(IncludeHandler& handler, std::ostream& stream,
		const std::string& indent = "\t", unsigned indentStartLevel = 0);
	virtual ~ConfigToTextVisitor();

	void setClassWritten(const ClassInfo* classInfo);
	bool isClassWritten(const ClassInfo* classInfo) const;

protected:
	void writeIndent();
};

}

}

#endif // CCfg_Io_ConfigToTextVisitor_h_INCLUDED
