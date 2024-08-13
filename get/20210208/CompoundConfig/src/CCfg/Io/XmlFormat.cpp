// =============================================================================
// XmlFormat.cpp
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

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include "XmlFormat.h"
#include "XmlToConfigVisitor.h"
#include "ConfigToXmlVisitor.h"
#include "StreamStorage.h"
#include "../Object.h"
#include "tinyxml.h"
#include <sstream>

using std::endl;

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlFormat class methods body
///////////////////////////////////////////////////////////////////////////////

const char XmlFormat::INDEX_XML_ATTR[]		= "id";
const char XmlFormat::TYPE_XML_ATTR[]		= "type";
const char XmlFormat::INCLUDE_XML_ATTR[]	= "include";
const char XmlFormat::UNIT_XML_ATTR[]		= "unit";
const char XmlFormat::RANGE_XML_ATTR[]		= "range";
const char* const XmlFormat::TYPE_VALUES[VALUETYPE_COUNT] = {"int", "hex", "real", "bool", "string"};

/**
 * Constructor.
 * @param useXmlDecl Flag specifying whether to use XML declarations.
 */
XmlFormat::XmlFormat(bool useXmlDecl) : useXmlDecl(useXmlDecl)
{
}

/**
 * Destructor.
 */
XmlFormat::~XmlFormat()
{
}

/**
 * Parses the storage data to produce a configuration tree.
 * @param storage the storage where is the data.
 * @param handler a callback to react to include directives.
 * @return Returns the created configuration tree.
 * @throws Storage::OpenError if there was an error opening the storage.
 */
std::auto_ptr<CCfg::Object> XmlFormat::parse(Storage& storage, IncludeHandler& handler)
{
	StreamStorage *stream = dynamic_cast<StreamStorage*>(&storage);
	assert(stream != 0);

	if(not stream->isOpened())
	{
		stream->open(false);
	}

	std::string data;
	stream->readAll(data);

	TiXmlDocument doc;
	doc.Parse(data.c_str());
	if(doc.Error())
	{
		throw SyntaxError(doc.ErrorId(), doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
	}

	XmlToConfigVisitor visitor(handler);
	doc.Accept(&visitor);
	return std::auto_ptr<CCfg::Object>(visitor.getBuiltConfig());
}

/**
 * Writes a configuration tree in a storage.
 * @param config the configuration tree.
 * @param storage the storage where data is written.
 * @param handler a callback to react to include directives.
 * @throws Storage::OpenError if there was an error opening the storage.
 */
void XmlFormat::format(const Object *config, Storage& storage, IncludeHandler& handler)
{
	StreamStorage *stream = dynamic_cast<StreamStorage*>(&storage);
	assert(stream != 0);

	if(not stream->isOpened())
	{
		stream->open(true);
	}

	std::ostringstream buffer;
	if (useXmlDecl)
	{
		buffer << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	}
	ConfigToXmlVisitor visitor(handler, buffer);
	config->accept(visitor);

	stream->write(buffer.str());
}

const char* XmlFormat::typeToString(ValueTypeId type)
{
	assert(type >= 0 and type < VALUETYPE_COUNT);
	return TYPE_VALUES[type];
}

ValueTypeId XmlFormat::typeFromString(const char* str)
{
	for(int i = 0; i < VALUETYPE_COUNT; ++i)
	{
		if(not strcmp(str, TYPE_VALUES[i]))
		{
			return static_cast<ValueTypeId>(i);
		}
	}
	return VALUETYPE_INVALID;
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlFormat::SyntaxError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param errorId the error id from the Xml parser.
 * @param errorDesc the error description from the Xml parser.
 * @param errorRow the row where the error occured in the Xml file.
 * @param errorCol the column where the error occured in the Xml file.
 */
XmlFormat::SyntaxError::SyntaxError(int errorId,
	const std::string& errorDesc,
	int errorRow, int errorCol)
	: errorId(errorId), errorDesc(errorDesc), errorRow(errorRow), errorCol(errorCol)
{
}

/**
 * Destructor.
 */
XmlFormat::SyntaxError::~SyntaxError() throw()
{
}

/**
 * Getter of the "errorId" property.
 * @return Returns the error id from the Xml parser.
 */
int XmlFormat::SyntaxError::getErrorId() const
{
	return errorId;
}

/**
 * Getter of the "errorDesc" property.
 * @return Returns the error description from the Xml parser.
 */
const std::string& XmlFormat::SyntaxError::getErrorDesc() const
{
	return errorDesc;
}

/**
 * Getter of the "errorRow" property.
 * @return Returns the row where the error occured in the Xml file.
 */
int XmlFormat::SyntaxError::getErrorRow() const
{
	return errorRow;
}

/**
 * Getter of the "errorCol" property.
 * @return Returns the column where the error occured in the Xml file.
 */
int XmlFormat::SyntaxError::getErrorCol() const
{
	return errorCol;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string XmlFormat::SyntaxError::formatMessage() const
{
	std::ostringstream message;
	message << "The text stream is not valid XML." << endl
			<< "At line " << errorRow << " column " << errorCol << ": "
			<< errorDesc << " (error #" << errorId << ").";
	return message.str();
}


///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlFormat::SemanticError class methods body
///////////////////////////////////////////////////////////////////////////////

const char* XmlFormat::SemanticError::errDescTable[ERROR_COUNT] =
	{"It cannot have both text and children XML elements.", // MIXING_TEXT_AND_ELEMENTS
	"The text of a scalar must not be split.",	// SPLIT_TEXT
	"It has an XML attribute that is not allowed by configuration specifications.", // INVALID_ATTRIBUTE
	"The scalar type is invalid."
	};

/**
 * Full constructor.
 * @param errorId the error id.
 * @param element the dump of the element where error occured.
 * @param errorRow the row where the error occured in the Xml file.
 * @param errorCol the column where the error occured in the Xml file.
 */
XmlFormat::SemanticError::SemanticError(Error errorId,
	const std::string& element, int errorRow, int errorCol)
	: errorId(errorId), element(element), errorRow(errorRow), errorCol(errorCol)
{
}

/**
 * Destructor.
 */
XmlFormat::SemanticError::~SemanticError() throw()
{
}

/**
 * Getter of the "errorId" property.
 * @return Returns the error id.
 */
int XmlFormat::SemanticError::getErrorId() const
{
	return errorId;
}

/**
 * Getter of the "element" property.
 * @return Returns the dump of the element where error occured.
 */
const std::string& XmlFormat::SemanticError::getElement() const
{
	return element;
}

/**
 * Getter of the "errorRow" property.
 * @return Returns the row where the error occured in the Xml file.
 */
int XmlFormat::SemanticError::getErrorRow() const
{
	return errorRow;
}

/**
 * Getter of the "errorCol" property.
 * @return Returns the column where the error occured in the Xml file.
 */
int XmlFormat::SemanticError::getErrorCol() const
{
	return errorCol;
}

/**
 * Retrieve the error description for a given id.
 * @return Returns the description of the error.
 */
const char* XmlFormat::SemanticError::getErrorDesc(Error id)
{
	unsigned index = static_cast<unsigned>(id);
	assert(index < ERROR_COUNT);
	return errDescTable[index];
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string XmlFormat::SemanticError::formatMessage() const
{
	std::ostringstream message;
	message << "The Xml document is not a valid configuration." << endl
			<< "The element: " << element
			<< " at line " << errorRow << " column " << errorCol
			<< " is causing a problem." << endl
			<< getErrorDesc(errorId);
	return message.str();
}

}

}
