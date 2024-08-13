// =============================================================================
// XmlToConfigVisitor.cpp
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

#include "XmlToConfigVisitor.h"
#include "IncludeHandler.h"
#include "../Attribute.h"
#include "../Object.h"
#include "../Scalar.h"
#include "../Range.h"
#include <cstring>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlToConfigVisitor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * Creates a visitor to read a configuration while ignoring included resources.
 */
XmlToConfigVisitor::XmlToConfigVisitor()
	: includeHandler(NoOpIncludeHandler::getInstance()),
	config(), currentAttr(0), pendingComment(), currentType(VALUETYPE_INVALID)
{
}

/**
 * Full constructor.
 * @param handler the callback to process includes.
 */
XmlToConfigVisitor::XmlToConfigVisitor(IncludeHandler& handler)
	: includeHandler(handler), config(), currentAttr(0), pendingComment()
{
}

/**
 * Destructor.
 */
XmlToConfigVisitor::~XmlToConfigVisitor()
{
}

/**
 * Getter of the "config" property.
 * @return Returns the configuration built while visiting the Xml Tree.
 */
std::auto_ptr<Object> XmlToConfigVisitor::getBuiltConfig()
{
	return config;
}

/**
 * Tests if an XML element stands for a configuration object or otherwise a
 * scalar.
 * @param element the XML element to analyse.
 * @return Returns true if it represents an object, false if it represents a
 * scalar.
 * @throws XmlFormat::SemanticError the element is neither a valid object nor a
 * valid scalar.
 */
bool XmlToConfigVisitor::isObject(const TiXmlElement& element)
{
	if(element.Attribute(XmlFormat::TYPE_XML_ATTR) != 0)
	{
		return false;
	}
	unsigned textChildrenCount = 0;
	unsigned elementChildrenCount = 0;

	const TiXmlNode* child = element.FirstChild();
	while(child != 0)
	{
		if(child->ToText())		textChildrenCount++;
		if(child->ToElement())	elementChildrenCount++;
		child = child->NextSibling();
	}

	if(textChildrenCount == 1 and elementChildrenCount == 0)
	{
		return false;
	}
	else if(textChildrenCount == 0)
	{
		return true;
	}
	else // invalid content
	{
		// preparing error data
		std::string elt = toString(element);
		int row = element.Row();
		int col = element.Column();

		if(textChildrenCount > 0 and elementChildrenCount > 0)
		{
			throw XmlSemError(XmlSemError::MIXING_TEXT_AND_ELEMENTS, elt, row, col);
		}
		else if(textChildrenCount > 1 and elementChildrenCount == 0)
		{
			throw XmlSemError(XmlSemError::SPLIT_TEXT, elt, row, col);
		}
		else // this is impossible but will prevent the compiler from yelling.
		{
			assert(false);
			return false;
		}
	}
}

/**
 * Extracts configuration object properties from an XML element.
 * @param element the XML element.
 * @param index the extracted index (out parameter).
 * @param include the extracted include path (out parameter).
 * @throws XmlFormat::SemanticError the element contains unrecognized
 * attributes.
 */
void XmlToConfigVisitor::extractObjectProperties(const TiXmlElement& element,
	std::string& index, std::string& include)
{
	const TiXmlAttribute* attr = element.FirstAttribute();
	while(attr != 0)
	{
		if(not strcmp(attr->Name(), XmlFormat::INDEX_XML_ATTR))
		{
			index = attr->Value();
		}
		else if(not strcmp(attr->Name(), XmlFormat::INCLUDE_XML_ATTR))
		{
			include = attr->Value();
		}
		else
		{
			throw XmlSemError(XmlSemError::INVALID_ATTRIBUTE, toString(element),
				element.Row(), element.Column());
		}
		attr = attr->Next();
	}
}

/**
 * Extracts configuration scalar properties from an XML element.
 * @param element the XML element.
 * @param [out] index  the extracted index.
 * @param [out] type the extracted type.
 * @param [out] unit the extracted unit.
 * @param [out] range the extracted range.
 * @throws XmlFormat::SemanticError the element contains unrecognized
 * attributes.
 */
void XmlToConfigVisitor::extractScalarProperties(const TiXmlElement& element,
	std::string& index, std::string& type, std::string& unit, std::string& range)
{
	const TiXmlAttribute* attr = element.FirstAttribute();
	while(attr != 0)
	{
		if(not strcmp(attr->Name(), XmlFormat::INDEX_XML_ATTR))
		{
			index = attr->Value();
		}
		else if(not strcmp(attr->Name(), XmlFormat::TYPE_XML_ATTR))
		{
			type = attr->Value();
		}
		else if(not strcmp(attr->Name(), XmlFormat::UNIT_XML_ATTR))
		{
			unit = attr->Value();
		}
		else if(not strcmp(attr->Name(), XmlFormat::RANGE_XML_ATTR))
		{
			range = attr->Value();
		}
		else
		{
			throw XmlSemError(XmlSemError::INVALID_ATTRIBUTE, toString(element),
				element.Row(), element.Column());
		}
		attr = attr->Next();
	}
}

/**
 * Converts an element opening markup to a string.
 * @param element the element to convert to string.
 * @return Returns the string representation of the element opening markup.
 */
std::string XmlToConfigVisitor::toString(const TiXmlElement& element) const
{
	struct NonRecursivePrinter : public TiXmlPrinter
	{
		NonRecursivePrinter() {}
		~NonRecursivePrinter() {}

		bool VisitEnter(const TiXmlElement& element, const TiXmlAttribute *firstAttribute)
		{
			TiXmlPrinter::VisitEnter(element, firstAttribute);
			return false;	// do not print content
		}
		bool VisitExit(const TiXmlElement& /*element*/)
		{
			// overridden method not called: do not print closing tag
			return false;	// do not print sibling elements
		}
		bool VisitEnter(const TiXmlDocument& /*doc*/) { return true; }
		bool VisitExit(const TiXmlDocument& /*doc*/) { return true; }
	};

	NonRecursivePrinter printer;
	printer.SetLineBreak("");
	element.Accept(&printer);
	return printer.Str();
}

/**
 * Callback reacting on the beginning of an XML document.
 * @param doc the XML document.
 * @return Returns always true ("go on recursing").
 */
bool XmlToConfigVisitor::VisitEnter(const TiXmlDocument& /*doc*/)
{
	return true;
}

/**
 * Callback reacting on the end of an XML document.
 * @param doc the XML document.
 * @return Returns always true ("go on recursing").
 */
bool XmlToConfigVisitor::VisitExit(const TiXmlDocument& /*doc*/)
{
	return true;
}

/**
 * Callback reacting on the beginning of an XML element.
 * @param element the XML element.
 * @param firstAttribute the first attribute of the element.
 * @return Returns always true ("go on recursing").
 * @throws XmlFormat::SemanticError the element does not respect
 * configuration format.
 * @throws Format::InvalidFormat an included configuration has syntaxic or
 * semantic errors.
 * @throws Storage::OpenError an included storage cannot be opened.
 */
bool XmlToConfigVisitor::VisitEnter(const TiXmlElement& element, const TiXmlAttribute* /*firstAttribute*/)
{
	if(isObject(element))
	{
		Object *newObj = new Object(element.Value());
		if(config.get() == 0)
		{
			config.reset(newObj);	// takes pointer ownership
		}
		else
		{
			assert(currentAttr != 0);
			currentAttr->append(newObj);	// takes pointer ownership
		}
		currentAttr = newObj;

		std::string index, include;
		extractObjectProperties(element, index, include);
		newObj->setIndex(index);
		newObj->setIncludeUrl(include);

		if(not include.empty())
		{
			includeHandler.parseInclude(*newObj);
		}
	}
	else
	{
		Scalar *newScalar = new Scalar(element.Value());
		currentAttr->append(newScalar);		// takes pointer ownership
		currentAttr = newScalar;

		std::string index, type, unit, range;
		extractScalarProperties(element, index, type, unit, range);
		newScalar->setIndex(index);
		if(not type.empty())
		{
			currentType = XmlFormat::typeFromString(type.c_str());
			if(currentType == VALUETYPE_INVALID)
			{
				throw XmlSemError(XmlSemError::INVALID_TYPE, toString(element),
					element.Row(), element.Column());
			}
		}
		else
		{
			currentType = VALUETYPE_INVALID; // set type as unspecified
		}

		if(not unit.empty())
		{
			newScalar->setUnit(unit);
		}
		if(not range.empty())
		{
			newScalar->setRange(Range::createFromString(range));
		}
	}
	if(not pendingComment.empty())
	{
		currentAttr->setPreComment(Comment(pendingComment));
		pendingComment.clear();
	}
	return true;
}

/**
 * Callback reacting on the end of an XML element.
 * @param element the XML element.
 * @return Returns always true ("go on recursing").
 */
bool XmlToConfigVisitor::VisitExit(const TiXmlElement& /*element*/)
{
	if(not pendingComment.empty())
	{
		currentAttr->setPostComment(Comment(pendingComment));
		pendingComment.clear();
	}
	currentAttr = currentAttr->getParent();
	return true;
}

/**
 * Callback reacting on the beginning of an XML text section.
 * @param text the XML text.
 * @return Returns always true ("go on recursing").
 */
bool XmlToConfigVisitor::Visit(const TiXmlText& text)
{
	if(currentType == VALUETYPE_INVALID) // if type is unspecified
	{
		currentAttr->setValue(text.Value()); // use type auto detection.
	}
	else
	{
		currentAttr->setValue(currentType, text.Value());
	}
	return true;
}

/**
 * Callback reacting on the beginning of an XML comment.
 * @param comment the XML comment.
 * @return Returns always true ("go on recursing").
 */
bool XmlToConfigVisitor::Visit(const TiXmlComment& comment)
{
	pendingComment += comment.Value();
	return true;
}

bool XmlToConfigVisitor::Visit(const TiXmlDeclaration& /*declaration*/)
{
	return true;
}

bool XmlToConfigVisitor::Visit(const TiXmlUnknown& /*unknown*/)
{
	return true;
}

}

}
