// =============================================================================
// ConfigToXmlVisitor.cpp
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

#include "ConfigToXmlVisitor.h"
#include "XmlFormat.h"
#include "IncludeHandler.h"
#include "../Scalar.h"
#include "../Object.h"
#include "../ScalarClassInfo.h"
#include "../Range.h"
#include "../Comment.h"
#include <cctype>

using std::endl;

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::ConfigToXmlVisitor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 * This constructor does not need an include handler.
 * It implies that all includes are ignored (they are not even written).
 */
ConfigToXmlVisitor::ConfigToXmlVisitor(std::ostream& stream,
	const std::string& indent, unsigned indentStartLevel)
	: ConfigToTextVisitor(NoOpIncludeHandler::getInstance(), stream, indent, indentStartLevel),
	emptyElement(false)
{
}

/**
 * Full constructor.
 * @param handler a callback to react to include directives.
 * @param stream the stream where to write XML.
 * @param indent the string defining one indent element (typicaly: n spaces or 1 tab)
 * @param indentStartLevel a number defining the indent level at which the root
 * element starts (no element will go below this indent level).
 */
ConfigToXmlVisitor::ConfigToXmlVisitor(
	IncludeHandler& handler, std::ostream& stream,
	const std::string& indent, unsigned indentStartLevel)
	: ConfigToTextVisitor(handler, stream, indent, indentStartLevel),
	emptyElement(false)
{
}

/**
 * Destructor.
 */
ConfigToXmlVisitor::~ConfigToXmlVisitor()
{
}

/**
 * Writes a comment in the XML stream.
 * Comment is indented and perform a line break at the end. The comment
 * content is written without rearangement. If this is a multiline comment
 * indentation won't be corrected.
 * @param comment the comment to write.
 */
void ConfigToXmlVisitor::writeComment(const Comment& comment)
{
	if(not comment.getText().empty())
	{
		writeIndent();
		stream << "<!--" << comment.getText() << "-->" << endl;
	}
}

/**
 * Callback function invoked when we enter a new object while visiting the
 * configuration.
 * @param objectCfg the object we enter.
 * @return Returns false to indicate that visitation should ignore nested
 * attributes, true otherwise. Typically, false is returned when an include
 * directive is processed, because attributes have already been saved
 * elsewhere.
 */
bool ConfigToXmlVisitor::onObjectEnter(const Object& objectCfg)
{
	writeComment(objectCfg.getPreComment());

	bool recurse = not objectCfg.empty();
	writeIndent();
	stream << '<' << objectCfg.getName();

	if(not objectCfg.getIndex().empty())
	{
		stream << ' ' << XmlFormat::INDEX_XML_ATTR
			<< "=\"" << XmlEscape(objectCfg.getIndex()) << '\"';
	}
	// ignore the include directive of the root object (because when formatting
	// an includer object, it would lead to infinite recursion).
	if(not objectCfg.getIncludeUrl().empty() and not nextIsRootObject)
	{
		// if the include was formatted in a separate storage.
		if(includeHandler.formatInclude(objectCfg))
		{
			// we have to put the include directive.
			stream << ' ' << XmlFormat::INCLUDE_XML_ATTR
				<< "=\"" << XmlEscape(objectCfg.getIncludeUrl()) << '\"';
			// there is no need to recurse.
			recurse = false;
		}
		else
		{
			// The include content will be inlined.
			// We should not put the include directive.
		}
	}
	if(recurse or not objectCfg.getPostComment().getText().empty())
	{
		stream << '>' << endl;
	}
	else
	{
		emptyElement = true;
	}

	nextIsRootObject = false;
	++indentLevel;

	return recurse;
}

/**
 * Callback function invoked when we leave an object while visiting the
 * configuration.
 * @param objectCfg the object we leave.
 * @return Always returns true (we never ignore siblings).
 */
bool ConfigToXmlVisitor::onObjectExit(const Object& objectCfg)
{
	--indentLevel;

	if(emptyElement)
	{
		stream << " />" << endl;
		emptyElement = false;
	}
	else
	{
		++indentLevel;
		writeComment(objectCfg.getPostComment());
		--indentLevel;
		writeIndent();
		stream << "</" << objectCfg.getName() << '>' << endl;
	}

	return true;
}

/**
 * Callback function invoked each time a scalar is encountered while
 * visiting the configuration.
 * @param scalar the scalar found.
 * @return Always returns true (we never ignore siblings).
 */
bool ConfigToXmlVisitor::onScalar(const Scalar& scalar)
{
	writeComment(scalar.getPreComment());

	writeIndent();
	stream << '<' << scalar.getName();

	// writes the xml attribute for the index
	if(not scalar.getIndex().empty())
	{
		stream << ' ' << XmlFormat::INDEX_XML_ATTR
			<< "=\"" << XmlEscape(scalar.getIndex()) << '\"';
	}
	// checks if the scalar value is ambiguous, and writes the xml attribute for the type if so
	Value parsed;
	const std::string& value = scalar.getValueAsString();
	parsed.setValue(value);
	if(value.empty() or parsed.getType() != scalar.getValueTypeId())
	{
		stream << ' ' << XmlFormat::TYPE_XML_ATTR
			<< "=\"" << XmlFormat::typeToString(scalar.getValueTypeId()) << '\"';
	}
	// writes class informations: unit and range if they exist.
	if(scalar.getClassInfo().get() != 0)
	{
		const ScalarClassInfo* classInfo = static_cast<ScalarClassInfo*>(scalar.getClassInfo().get());
		if(not isClassWritten(classInfo))
		{
			if(not classInfo->getUnit().empty())
			{
				stream << ' ' << XmlFormat::UNIT_XML_ATTR
					<< "=\"" << XmlEscape(classInfo->getUnit()) << '\"';
			}
			if(classInfo->getRange() != 0)
			{
				stream << ' ' << XmlFormat::RANGE_XML_ATTR
					<< "=\"" << XmlEscape(classInfo->getRange()->toString()) << '\"';
			}
			setClassWritten(classInfo);
		}
	}
	// writes the scalar value.
	stream << '>' << XmlEscape(value);
	// writes the post-comment.
	if(not scalar.getPostComment().getText().empty())
	{
		++indentLevel;
		stream << endl;
		writeComment(scalar.getPostComment());
		--indentLevel;
		writeIndent();
	}
	stream << "</" << scalar.getName() << '>' << endl;

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::ConfigToXmlVisitor::XmlEscape class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 *
 */
std::ostream& operator<<(std::ostream& stream, const ConfigToXmlVisitor::XmlEscape& esc)
{
	for(size_t i = 0; i < esc.str.length(); ++i)
	{
		char c = esc.str[i];
		if(ispunct(c))
		{
			switch(c)
			{
			case '&':stream << "&amp;";
				break;
			case '<':stream << "&lt;";
				break;
			case '>':stream << "&gt;";
				break;
			case '"':stream << "&quot;";
				break;
			case '\'':stream << "&apos;";
				break;
			default:
				stream << c;
			}
		}
		else if(isprint(c) or c == '\t' or c == '\n' or c == '\r')
		{
			stream << c;
		}
		else
		{
			stream << "&#x" << std::hex << static_cast<int>(c) << std::dec << ';';
		}
	}
	return stream;
}

}
}
