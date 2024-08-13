// =============================================================================
// XmlToConfigVisitor.h
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

#ifndef CCfg_Io_XmlToConfigVisitor_h_INCLUDED
#define CCfg_Io_XmlToConfigVisitor_h_INCLUDED

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include "XmlFormat.h"
#include "tinyxml.h"
#include <string>
#include <memory>

namespace CCfg
{
// forward declarations
class Attribute;

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class builds a configuration while visiting a TinyXml tree.
 *
 * @author fchateau
 */
class XmlToConfigVisitor : public TiXmlVisitor
{
	IncludeHandler& includeHandler;	///< the callback to process includes.
	std::auto_ptr<Object> config;	///< the configuration built while visiting the Xml Tree.

	// transient attributes
	Attribute* currentAttr;	///< the config attribute for currently opened Xml element.
	/// a concatenation of comments waiting to be associated to a configuration attribute.
	std::string pendingComment;
	ValueTypeId currentType;

public:
	// construction/destruction
	XmlToConfigVisitor();
	explicit XmlToConfigVisitor(IncludeHandler& handler);
	virtual ~XmlToConfigVisitor();

	// getters
	std::auto_ptr<Object>  getBuiltConfig();

	// visitation
	bool VisitEnter(const TiXmlDocument& doc);
	bool VisitExit(const TiXmlDocument& doc);

	bool VisitEnter(const TiXmlElement& element, const TiXmlAttribute* firstAttribute);
	bool VisitExit(const TiXmlElement& element);

	bool Visit(const TiXmlText& text);
	bool Visit(const TiXmlComment& comment);

	bool Visit(const TiXmlDeclaration& declaration);
	bool Visit(const TiXmlUnknown& unknown);

private:
	bool isObject(const TiXmlElement& element);
	void extractObjectProperties(const TiXmlElement& element,
			std::string& index, std::string& include);
	void extractScalarProperties(const TiXmlElement& element,
			std::string& index, std::string& type,
			std::string& unit, std::string& range);

	std::string toString(const TiXmlElement& element) const;
};

}

}

#endif // CCfg_Io_XmlToConfigVisitor_h_INCLUDED
