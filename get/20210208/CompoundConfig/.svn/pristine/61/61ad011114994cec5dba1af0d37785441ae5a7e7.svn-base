// =============================================================================
// XmlFormat.h
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

#ifndef CCfg_Io_XmlFormat_h_INCLUDED
#define CCfg_Io_XmlFormat_h_INCLUDED

#include "TextFormat.h"
#include "../TypeInfo.h"

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class provide methods to map the configuration semantics in XML.
 *
 * author fchateau
 */
class XmlFormat : public TextFormat
{
public:
	// inner classes
	class SyntaxError;
	class SemanticError;

	// construction/destruction
	XmlFormat(bool useXmlDecl=true);
	virtual ~XmlFormat();

	std::auto_ptr<CCfg::Object> parse(Storage& storage, IncludeHandler& handler);
	void format(const Object *config, Storage& storage, IncludeHandler& handler);

	static const char INDEX_XML_ATTR[];
	static const char TYPE_XML_ATTR[];
	static const char INCLUDE_XML_ATTR[];
	static const char UNIT_XML_ATTR[];
	static const char RANGE_XML_ATTR[];
	static const char* const TYPE_VALUES[VALUETYPE_COUNT];

	static const char* typeToString(ValueTypeId type);
	static ValueTypeId typeFromString(const char* str);
private:
	bool useXmlDecl; ///< Whether to start documents with an XML declaration.
};



/**
 * Reports that the XML stream does not respect Xml specifications.
 */
class XmlFormat::SyntaxError : public Format::SyntaxError
{
	int errorId;	///< the error id from the Xml parser.
	std::string errorDesc;	///< the error description from the Xml parser.
	int errorRow;	///< the row where the error occured in the Xml file.
	int errorCol;	///< the column where the error occured in the Xml file.

public:
	SyntaxError(int errorId, const std::string& errorDesc, int errorRow, int errorCol);
	~SyntaxError() throw();

	int getErrorId() const;
	const std::string& getErrorDesc() const;
	int getErrorRow() const;
	int getErrorCol() const;

private:
	std::string formatMessage() const;
};

/**
 * Reports that the XML tree does not respect the semantics of our
 * configuration.
 */
class XmlFormat::SemanticError : public Format::SemanticError
{
public:
	enum Error {MIXING_TEXT_AND_ELEMENTS, SPLIT_TEXT, INVALID_ATTRIBUTE, INVALID_TYPE, ERROR_COUNT};

	SemanticError(Error errorId, const std::string& element, int errorRow, int errorCol);
	~SemanticError() throw();

	int getErrorId() const;
	const std::string& getElement() const;
	int getErrorRow() const;
	int getErrorCol() const;

	static const char* getErrorDesc(Error id);

private:
	std::string formatMessage() const;

	Error errorId;	///< the error id.
	std::string element;	///< the dump of the element where error occured.
	int errorRow;	///< the row where the error occured in the Xml file.
	int errorCol;	///< the column where the error occured in the Xml file.

	static const char* errDescTable[ERROR_COUNT];
};

typedef XmlFormat::SemanticError XmlSemError;

}

}

#endif // CCfg_Io_XmlFormat_h_INCLUDED
