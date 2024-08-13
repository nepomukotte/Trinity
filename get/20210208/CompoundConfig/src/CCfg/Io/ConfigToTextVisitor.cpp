// =============================================================================
// ConfigToTextVisitor.cpp
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

#include "ConfigToTextVisitor.h"
#include <algorithm>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::ConfigToTextVisitor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param handler a callback to react to include directives.
 * @param stream the stream where to write text.
 * @param indent the string defining one indent element (typically: n spaces or 1 tab)
 * @param indentStartLevel a number defining the indent level at which the root
 * element starts (no element will go below this indent level).
 */
ConfigToTextVisitor::ConfigToTextVisitor(IncludeHandler& handler, std::ostream& stream,
	const std::string& indent, unsigned indentStartLevel)
	: includeHandler(handler), stream(stream),
	indentPattern(indent), writtenClasses(),
	indentLevel(indentStartLevel), nextIsRootObject(true)
{
}

/**
 * Destructor.
 */
ConfigToTextVisitor::~ConfigToTextVisitor()
{
}

/**
 * Writes the indentation according to its current level.
 */
void ConfigToTextVisitor::writeIndent()
{
	for(unsigned i = 0; i < indentLevel; i++)
	{
		stream << indentPattern;
	}
}

/**
 * Registers the fact that the specified ClassInfo has already been written to
 * the current file.
 */
void ConfigToTextVisitor::setClassWritten(const ClassInfo* classInfo)
{
	writtenClasses.push_back(classInfo);
}

/**
 * Asks whether the specified ClassInfo has already been written to the current file.
 */
bool ConfigToTextVisitor::isClassWritten(const ClassInfo* classInfo) const
{
	return find(writtenClasses.begin(), writtenClasses.end(), classInfo) != writtenClasses.end();
}

}
}
