// =============================================================================
// ConfigToXmlVisitor.h
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

#ifndef CCfg_Io_ConfigToXmlVisitor_h_INCLUDED
#define CCfg_Io_ConfigToXmlVisitor_h_INCLUDED

#include "ConfigToTextVisitor.h"
#include <string>
#include <ios>

namespace CCfg
{
// forward declarations
class Comment;

namespace Io
{
// forward declarations
class IncludeHandler;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class outputs a configuration tree to a XML stream.
 *
 * @author fchateau
 */
class ConfigToXmlVisitor : public ConfigToTextVisitor
{
	// transient attributes
	bool emptyElement;

public:
	// inner class
	struct XmlEscape;

	// construction/destruction
	explicit ConfigToXmlVisitor(std::ostream& stream,
		const std::string& indent = "\t", unsigned indentStartLevel = 0);
	ConfigToXmlVisitor(IncludeHandler& handler, std::ostream& stream,
		const std::string& indent = "\t", unsigned indentStartLevel = 0);
	virtual ~ConfigToXmlVisitor();

	// visitation
	virtual bool onObjectEnter(const Object& objectCfg);
	virtual bool onObjectExit(const Object& objectCfg);
	virtual bool onScalar(const Scalar& scalar);

protected:
	void writeComment(const Comment& comment);
};



/**
 * Stream manipulator to make a string compliant with Xml character restrictions.
 */
struct ConfigToXmlVisitor::XmlEscape
{
	const std::string& str;

	XmlEscape(const std::string& str);
};

std::ostream& operator<<(std::ostream& stream, const ConfigToXmlVisitor::XmlEscape& esc);






////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::ConfigToXmlVisitor::XmlEscape inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param str the string to escape before outputing to an XML stream.
 */
inline ConfigToXmlVisitor::XmlEscape::XmlEscape(const std::string& str)
	: str(str)
{
}

}
}

#endif // CCfg_Io_ConfigToXmlVisitor_h_INCLUDED
