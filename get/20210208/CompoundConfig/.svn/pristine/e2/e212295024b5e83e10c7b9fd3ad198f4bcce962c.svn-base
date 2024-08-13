// =============================================================================
// Attribute.cpp
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

#include "CConfig.h"
#include "View/InheritedIterator.h"
#include "View/InheritersIterator.h"
#include "View/ClassIterator.h"
#include "View/Search.h"
#include "Object.h"
#include "Range.h"
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

using std::vector;
using std::string;
using std::ostringstream;
using std::make_pair;

#undef max

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
CConfig::CConfig()
	: baseAttribute(0), path(), baseDepth(0), backtrackLevel(0)
{
}

/**
 * Transtyping constructor: initializes a config for the specified attribute.
 * @param baseAttribute the attribute defining the entry point to the
 * configuration and the initial path.
 * Must not be null.
 */
CConfig::CConfig(Attribute* baseAttribute)
	// a little hack to assert before something bad happends
	: baseAttribute( (assert(baseAttribute != 0), baseAttribute) ),
	path(baseAttribute->getPath()),	// we use the path of the base attribute.
	baseDepth(path.size()-1),	// the base attribute is at the top of the path.
	backtrackLevel(0)	// the path directly designate the base attribute.
{
}

/**
 * Constructor: initializes a config for the specified attribute, with a
 * different path.
 * \note
 * By default the path is the base attribute path, but in some cases taking
 * the baseAttribute path would lead to a loss of information.
 * Indeed, if baseAttribute is contained in one or more default objects at some
 * levels of its ascendancy, the index that was initially searched may be lost.
 * This constructor exist precisely to cope with this issue. It allows to
 * define a path with relevant indexes.
 * @param baseAttribute an attribute matching the path.
 * @param path the path that the virtual attribute represented by the config.
 * @param backtrackLevel The number of backtrackings when searching the base
 * attribute from the root, using the path.
 * should have.
 */
CConfig::CConfig(Attribute* baseAttribute, const Path& path,
	int backtrackLevel)
	: baseAttribute(baseAttribute),
	path(path),
	baseDepth(path.size()-1),
	backtrackLevel(backtrackLevel)
{
	assert(baseAttribute != 0);
}

/**
 * Copy constructor: initializes a config with the resolution of the source config.
 * If the source config has been altered with operator() or find, when we copy
 * it, the destination will be initialized to the resolution of the config.
 * @param config the config to copy.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::CConfig(const CConfig& config)
	: baseAttribute(),
	path(config.getPath()),
	baseDepth(path.size()-1),
	backtrackLevel(config.getBacktrackLevel())
{
	if(config.baseAttribute != 0)
	{
		// resolve must be called after copying the path, because it calls resetPath.
		View::Result res = config.resolve();
		baseAttribute = res.getAttribute();
		backtrackLevel = std::max(backtrackLevel, res.getBacktrackLevel());
	}
}

/**
 * Destructor.
 */
CConfig::~CConfig()
{
}

/**
 * Assignment operator: overwrites a config with the resolution of another.
 * @see CCfg::CConfig::CConfig(const CConfig&)
 * @param config the config to copy.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig& CConfig::operator =(const CConfig& config)
{
	path = config.getPath();
	baseDepth = path.size()-1;
	if(config.baseAttribute != 0)
	{
		// resolve must be called after copying the path, because it calls
		// resetPath.
		View::Result res = config.resolve();
		baseAttribute = res.getAttribute();
		backtrackLevel = std::max(res.getBacktrackLevel(), config.getBacktrackLevel());
	}
	return *this;
}

/**
 * Setter of the property "baseAttribute".
 * Resets the config to make it point on the specified attribute.
 * @param newBase the attribute on which the config will point from now on.
 */
void CConfig::setBaseAttribute(Attribute* newBase)
{
	assert(newBase != 0);
	baseAttribute = newBase;
	path = baseAttribute->getPath();
	baseDepth = path.size()-1;
	backtrackLevel = 0;
}

/**
 * Creates each object or scalar of the config path that doesn't exist in the model.
 * @param lastInstType indicates whether the tail of the path represents an
 * scalar or an object.
 * All other path nodes will be created as object (which is logical, a scalar
 * would mark the end of the path).
 * @throw CCfg::UnsupportedOperation if a path element exists but is a scalar
 * instead of an object.
 */
Attribute& CConfig::create(InstanceType lastInstType)
{
	Attribute* root = baseAttribute->getRootConfig();
	Path pathFromRoot = path.subPath(1);

	Attribute& newAttribute = root->require(lastInstType, pathFromRoot);
	resetPath();
	return newAttribute;
}

/**
 * Returns the CCfg::Attribute designated by the current config path.
 * It may return 0 if no one matches the path.
 */
Attribute* CConfig::search() const
{
	assert(baseAttribute != 0);

	View::Search<Attribute> s;
	Attribute* attribute = s.search(path, *baseAttribute, baseDepth);
	resetPath();
	return attribute;
}

/**
 * Returns the CCfg::Attribute designated by the current config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
Attribute& CConfig::find() const
{
	assert(baseAttribute != 0);

	View::Search<Attribute> s;
	Attribute* attribute = s.search(path, *baseAttribute, baseDepth);
	if(attribute == 0)
	{
		Path invalidPath = path;
		resetPath();
		throw NotFound(invalidPath);
	}
	resetPath();
	return *attribute;
}

/**
 * Computes the attribute matching the config path, starting from the
 * baseAttribute.
 * @return Returns the resolution result: the attribute matching the config and
 * the number of backtracks.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
View::Result CConfig::resolve() const
{
	assert(baseAttribute != 0);

	View::Search<Attribute> s;
	Attribute* attribute = s.search(path, *baseAttribute, baseDepth);
	if(attribute == 0)
	{
		throw NotFound(path);
	}
	resetPath();
	return View::Result(attribute, std::max(backtrackLevel, s.getMaxBacktrackLevel()));
}

/**
 * Alias for the append method (const char version).
 */
CConfig& CConfig::operator ()(const char* name, const char* index)
{
	return append(name, index);
}

/**
 * Alias for the append method (const char version).
 */
CConfig& CConfig::operator ()(const char* name, int index)
{
	return append(name, index);
}

/**
 * Alias for the appendPath method.
 * @see appendPath(const std::string&)
 */
CConfig& CConfig::at(const char* relativePath)
{
	return appendPath(relativePath);
}

/**
 * Temporarily appends some attribute identifiers to the current path of the config.
 * The current path is reset to its initialization value each time a method
 * other than a path related one is called.
 */
CConfig& CConfig::append(const std::string& name, const std::string& index)
{
	path << PathNode(name, index);
	return *this;
}

/**
 * Temporarily appends some attribute identifiers to the current path of the config.
 * The current path is reset to its initialization value each time a method
 * other than a path related one is called.
 */
CConfig& CConfig::append(const std::string& name, int index)
{
	path << PathNode(name, index);
	return *this;
}

/**
 * Temporarily appends a path to the current path of the config.
 * The current path is reset to its initialization value each time a method
 * other than a path related one is called.
 */
CConfig& CConfig::appendPath(const Path& relativePath)
{
	path.append(relativePath);
	return *this;
}

/**
 * Temporarily appends a path to the current path of the config.
 * The current path is reset to its initialization value each time a method
 * other than a path related one is called.
 */
CConfig& CConfig::appendPath(const std::string& relativePath)
{
	return appendPath(Path(relativePath));
}

/**
 * Creates an iterator over the attributes of the object pointed by the config.
 * It is the result of the merging of every attributes of every default
 * object of the object pointed by the config.
 * @param filter enable to ignore attribute whose name does not match the filter.
 * For the moment, only exact name matching is supported.
 * @return Returns an iterator over the object's attributes.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::Iterator CConfig::iterate(std::string filter, bool includeDefault) const
{
	return CConfig::Iterator(*this, filter, includeDefault);
}

/**
 * Creates an iterator over the default objects of the object pointed by the config.
 * \note
 * It works also for attributes.
 * @return Returns a default objects iterator.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
View::InheritedIterator CConfig::getInheritedAttributes() const
{
	Path parentPath = path; // we must backup the path, because resolve will reset it.
	return View::InheritedIterator(resolve(), parentPath);
}

/**
 * Creates an iterator over the inheriting objects of the object pointed by the config.
 * @return Returns a heir objects iterator.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
View::InheritersIterator CConfig::getInheriterAttributes() const
{
	Path parentPath = path; // we must backup the path, because resolve will reset it.
	return View::InheritersIterator(resolve(), parentPath);
}

View::ClassIterator CConfig::getClassInstances() const
{
	Path attrPath = path;
	resetPath();
	return View::ClassIterator(*baseAttribute->getRootConfig(), attrPath);
}

/**
 * Resets the path to its initial value, before being altered by searches.
 */
void CConfig::resetPath() const
{
	while(path.size()-1 > baseDepth)
	{
		path.removeTail();
	}
}

/**
 * Returns the name of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
const std::string& CConfig::getName() const
{
	return find().getName();
}

/**
 * Returns the index of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
const std::string& CConfig::getIndex() const
{
	return find().getIndex();
}

/**
 * Returns the pre-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
Comment& CConfig::getPreComment()
{
	return find().getPreComment();
}

/**
 * Returns the post-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
Comment& CConfig::getPostComment()
{
	return find().getPostComment();
}

/**
 * Returns the pre-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
const Comment& CConfig::getPreComment() const
{
	return find().getPreComment();
}

/**
 * Returns the post-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
const Comment& CConfig::getPostComment() const
{
	return find().getPostComment();
}

/**
 * Returns the include URL of the configuration object under this config.
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws UnsupportedOperation if the configuration is not a CCfg::Object.
 */
const std::string& CConfig::getIncludeUrl() const
{
	return find().getIncludeUrl();
}

/**
 * Changes the name of the current configuration attribute.
 */
void CConfig::setName(const std::string& value)
{
	if(hasBasePath())
	{
		baseAttribute->setName(value);
		path.getTail().name = value;
	}
	else
	{
		find().setName(value);
	}
}

/**
 * Changes the index of the current configuration attribute.
 */
void CConfig::setIndex(const std::string& value)
{
	if(hasBasePath())
	{
		baseAttribute->setIndex(value);
		path.getTail().index = value;
	}
	else
	{
		find().setIndex(value);
	}
}

/**
 * Defines the pre-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
void CConfig::setPreComment(const Comment& value)
{
	find().setPreComment(value);
}

/**
 * Defines the post-comment of the attribute referenced by the config path.
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
void CConfig::setPostComment(const Comment& value)
{
	find().setPostComment(value);
}

/**
 * Materializes the config and defines the include URL of the leaf object.
 * @throw CCfg::UnsupportedOperation if a path element exists but is a scalar
 * instead of an object.
 */
void CConfig::setIncludeUrl(const std::string& value)
{
	create(OBJECT_INST).setIncludeUrl(value);
}

/**
 * Assigns a value to the scalar designated by the config path.
 * \note
 * If any object or attribute in the path does not exist, we create it.
 * It means that assignation modify the model structure instead of resolving
 * the config and altering the matching scalar.
 * @param value the value to put in the scalar.
 * @return Returns the assigned attribute.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 */
template<typename T> Attribute& CConfig::assign(T value)
{
	Attribute& attr = create(SCALAR_INST);
	attr = value;
	return attr;
}

/**
 * Assigns a string to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(const String&)
 */
Attribute& CConfig::operator =(const String& value)
{
	return assign(value);
}

/**
 * Assigns an integer to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Int)
 */
Attribute& CConfig::operator =(Int value)
{
	return assign(value);
}

/**
 * Assigns an hexadecimal to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Hex)
 */
Attribute& CConfig::operator =(Hex value)
{
	return assign(value);
}

/**
 * Assigns a real number to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Real)
 */
Attribute& CConfig::operator =(Real value)
{
	return assign(value);
}

/**
 * Assigns a boolean to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Bool)
 */
Attribute& CConfig::operator =(Bool value)
{
	return assign(value);
}

/**
 * Assigns a string to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(const char*)
 */
Attribute& CConfig::operator =(const char* value)
{
	return assign(value);
}

/**
 * Assigns a generic value to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(const Value&)
 */
Attribute& CConfig::operator =(const Value& value)
{
	return assign(value);
}

/**
 * Assigns an integer to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Int)
 */
Attribute& CConfig::operator =(short value)
{
	return assign(value);
}

/**
 * Assigns an integer to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Int)
 */
Attribute& CConfig::operator =(long value)
{
	return assign(value);
}

/**
 * Assigns an hexadecimal to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Hex)
 */
Attribute& CConfig::operator =(unsigned short value)
{
	return assign(value);
}

/**
 * Assigns an hexadecimal to the concrete attribute.
 * @see CCfg::CConfig::assign(T)
 * @see Attribute::operator =(Hex)
 */
Attribute& CConfig::operator =(unsigned long value)
{
	return assign(value);
}

/**
 * Returns the string value in the concrete attribute.
 * @see Attribute::operator const String&() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator const String&() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the integer value in the concrete attribute.
 * @see Attribute::operator Int() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator Int() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the hexadecimal value in the concrete attribute.
 * @see Attribute::operator Hex() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator Hex() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the real number value in the concrete attribute.
 * @see Attribute::operator Real() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator Real() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the boolean value in the concrete attribute.
 * @see Attribute::operator Bool() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator Bool() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the string value in the concrete attribute.
 * @see Attribute::operator const char*() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator const char*() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the integer value in the concrete attribute.
 * @see Attribute::operator Int() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator short() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the integer value in the concrete attribute.
 * @see Attribute::operator Int() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator long() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the hexadecimal value in the concrete attribute.
 * @see Attribute::operator Hex() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator unsigned short() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the hexadecimal value in the concrete attribute.
 * @see Attribute::operator Hex() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 */
CConfig::operator unsigned long() const
{
	Attribute& attr = find();
	return attr;
}

/**
 * Returns the unit of the concrete attribute.
 * @see Attribute::getUnit() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
const std::string& CConfig::getUnit() const
{
	return find().getUnit();
}

/**
 * Returns the range of a Scalar.
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Range* CConfig::getRange() const
{
	return find().getRange();
}

/**
 * Materializes the scalar designated by the config path, and set its unit.
 * @see Attribute::setUnit(const std::string&)
 * @see Attribute& create(InstanceType)
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 */
void CConfig::setUnit(const std::string& value)
{
	create(SCALAR_INST).setUnit(value);
}

/**
 * Defines the range of the equivalence class of this Scalar.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 */
void CConfig::setRange(RangePtr range)
{
	create(SCALAR_INST).setRange(range);
}

/**
 * Defines the range of the equivalence class of this Scalar.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 */
void CConfig::setRange(const Range& range)
{
	create(SCALAR_INST).setRange(range);
}

/**
 * Returns the value of the concrete attribute.
 * @see Attribute::getValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
const Value& CConfig::getValue() const
{
	return find().getValue();
}

/**
 * Returns the type id of the value in the concrete attribute.
 * @see Attribute::getValueTypeId() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
ValueTypeId CConfig::getValueTypeId() const
{
	return find().getValueTypeId();
}

/**
 * Returns the string value in the concrete attribute.
 * @see Attribute::getStringValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
const String& CConfig::getStringValue() const
{
	return find().getStringValue();
}

/**
 * Returns the integer value in the concrete attribute.
 * @see Attribute::getIntValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Int CConfig::getIntValue() const
{
	return find().getIntValue();
}

/**
 * Returns the hexadecimal value in the concrete attribute.
 * @see Attribute::getHexValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Hex CConfig::getHexValue() const
{
	return find().getHexValue();
}

/**
 * Returns the real number value in the concrete attribute.
 * @see Attribute::getRealValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Real CConfig::getRealValue() const
{
	return find().getRealValue();
}

/**
 * Returns the boolean value in the concrete attribute.
 * @see Attribute::getBoolValue() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Bool CConfig::getBoolValue() const
{
	return find().getBoolValue();
}

/**
 * Returns the string value in the concrete attribute if it exists. If not, returns the specified default value.
 * @see Attribute::getStringValue() const
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
const String& CConfig::getStringValue(const String& defaultValue) const
{
	try
	{
		return find().getStringValue();
	}
	catch (const CConfig::NotFound& e)
	{
		return defaultValue;
	}
}

/**
 * Returns the integer value in the concrete attribute if it exists. If not, returns the specified default value.
 * @see Attribute::getIntValue() const
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Int CConfig::getIntValue(Int defaultValue) const
{
	try
	{
		return find().getIntValue();
	}
	catch (const CConfig::NotFound& e)
	{
		return defaultValue;
	}
}

/**
 * Returns the hexadecimal value in the concrete attribute if it exists. If not, returns the specified default value.
 * @see Attribute::getHexValue() const
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Hex CConfig::getHexValue(Hex defaultValue) const
{
	try
	{
		return find().getHexValue();
	}
	catch (const CConfig::NotFound& e)
	{
		return defaultValue;
	}
}

/**
 * Returns the real value in the concrete attribute if it exists. If not, returns the specified default value.
 * @see Attribute::getRealValue() const
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Real CConfig::getRealValue(Real defaultValue) const
{
	try
	{
		return find().getRealValue();
	}
	catch (const CConfig::NotFound& e)
	{
		return defaultValue;
	}
}

/**
 * Returns the boolean value in the concrete attribute if it exists. If not, returns the specified default value.
 * @see Attribute::getBoolValue() const
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Bool CConfig::getBoolValue(Bool defaultValue) const
{
	try
	{
		return find().getBoolValue();
	}
	catch (const CConfig::NotFound& e)
	{
		return defaultValue;
	}
}

/**
 * Returns the conversion to string of the value in the concrete attribute.
 * @see Attribute::getValueAsString() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
const String& CConfig::getValueAsString() const
{
	return find().getValueAsString();
}

/**
 * Returns the conversion to integer of the value in the concrete attribute.
 * @see Attribute::getValueAsInt() const
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Int CConfig::getValueAsInt() const
{
	return find().getValueAsInt();
}

/**
 * Returns the conversion to hexadecimal of the value in the concrete attribute.
 * @see Attribute::getValueAsHex()
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Hex CConfig::getValueAsHex() const
{
	return find().getValueAsHex();
}

/**
 * Returns the conversion to real number of the value in the concrete attribute.
 * @see Attribute::getValueAsReal()
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Real CConfig::getValueAsReal() const
{
	return find().getValueAsReal();
}

/**
 * Returns the conversion to boolean of the value in the concrete attribute.
 * @see Attribute::getValueAsBool()
 * @throws CConfig::NotFound if no attribute matches the config path.
 * @throws CCfg::UnsupportedOperation if the concrete attribute is not a scalar.
 */
Bool CConfig::getValueAsBool() const
{
	return find().getValueAsBool();
}

/**
 * Materializes the scalar designated by the config path, and set its value.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setValue(const std::string&)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setValue(const std::string& value)
{
	create(SCALAR_INST).setValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setValue(const Value&)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setValue(const Value& value)
{
	create(SCALAR_INST).setValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as a
 * string.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setStringValue(const String&)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setStringValue(const String& value)
{
	create(SCALAR_INST).setStringValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as an
 * integer.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setIntValue(Int)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setIntValue(Int value)
{
	create(SCALAR_INST).setIntValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as an
 * hexadecimal.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setHexValue(Hex)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setHexValue(Hex value)
{
	create(SCALAR_INST).setHexValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as a
 * real number.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setRealValue(Real)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setRealValue(Real value)
{
	create(SCALAR_INST).setRealValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as a
 * boolean.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setBoolValue(Bool)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setBoolValue(Bool value)
{
	create(SCALAR_INST).setBoolValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as an
 * integer.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setIntValue(const String&)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setIntValue(const String& value)
{
	create(SCALAR_INST).setIntValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as an
 * hexadecimal.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setHexValue(Hex)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setHexValue(const String& value)
{
	create(SCALAR_INST).setHexValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as a
 * real number.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setRealValue(Real)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setRealValue(const String& value)
{
	create(SCALAR_INST).setRealValue(value);
}

/**
 * Materializes the scalar designated by the config path, and set its value as a
 * boolean.
 * @throw CCfg::UnsupportedOperation the model is incompatible with the config
 * path.
 * @see Attribute::setBoolValue(Bool)
 * @see Attribute& create(InstanceType)
 */
void CConfig::setBoolValue(const String& value)
{
	create(SCALAR_INST).setBoolValue(value);
}


/**
 * Make the list of attributes of this config.
 * @param nameFilter enable to ignore attribute whose name does not match the filter.
 * @return Returns the list of attributes of this config.
 */
ConfigList CConfig::getAttributes(const std::string& nameFilter, bool includeDefault) const
{
	Path basePath = getPath(); // back up the path before config resolution.
	ConfigList outViews;
	// First, we enumerate inherited objects of the object designated by the
	// config.
	InheritedList defaults;
	listInherited(defaults);

	// Then, we enumerate every attribute of every inherited object.
	AttributeList attrNames;
	AttributeLookup attributeLookup;
	listAttributes(defaults, nameFilter, attrNames, attributeLookup);

	// Finally, we create a config for each attribute found.
	// This pass cannot be grouped with the previous, because:
	// 1/ we need attributes to be ordered
	// 2/ the config should point to the nearest overriden attribute
	for(AttributeList::iterator iAttr = attrNames.begin(),
		endAttr = attrNames.end(); iAttr != endAttr; ++iAttr)
	{
		// get the search result for the current attribute name.
		View::Result& res = attributeLookup[*iAttr];

		// construct the path of the new config using the path of the parent
		// config as stem. (we cannot take the path of the attribute itself,
		// because we would lose information).
		Path attrPath = basePath;
		attrPath.append(res);
		// creates a config and add it to the list
		if (includeDefault or res->getIndex() != CCfg::DEFAULT_OBJ_INDEX)
		{
			outViews.push_back(CConfig(res, attrPath, res.getBacktrackLevel()));
		}
	}
	return outViews;
}

/**
 * Makes the list of objects matching the path (the parent and its default
 * objects).
 * @param defaults [out] a vector that will be filled with pointers to every
 * objects related to the parent object by inheritance (including the parent
 * itself).
 */
void CConfig::listInherited(InheritedList& defaults) const
{
	CConfig resolvedAttr = *this;
	View::InheritedIterator def = resolvedAttr.getInheritedAttributes();

	defaults.push_back(resolvedAttr.resolve());
	while(def.hasNext())
	{
		def.next();
		defaults.push_back(def.get());
	}
}

/**
 * Makes the list of every attributes (removing doublets) for each object in
 * the parent inherited objects list.
 * @param defaults [in] the list of inherited object of the parent (including itself).
 * @param filter [in] enable to ignore attribute whose name does not match the filter.
 * @param attributes [out] the list of attributes found, in the order
 * they appear in the configuration.
 * @param attributeLookup [out] a map that associate the node name of each attribute
 * to a pointer to the parent's nearest object (with respect to configuration
 * priorities).
 */
void CConfig::listAttributes(const InheritedList& defaults,
	const std::string& filter,
	AttributeList& attributes,
	AttributeLookup& attributeLookup)
{
	// Beginning by the last default object, we iterate on these objects.
	// This order is used to preserve the order in which attributes are
	// declared in the configuration file.
	for(InheritedList::const_reverse_iterator iDefault = defaults.rbegin(),
		endDefault = defaults.rend(); iDefault != endDefault; ++iDefault)
	{
		// for each attribute of the default object.
		for(Attribute::const_iterator iAttr = (*iDefault)->begin(),
			endAttr = (*iDefault)->end(); iAttr != endAttr; ++iAttr)
		{
			// only add attribute when there is no filter
			// or when its name matched the filter
			if(filter.empty() or iAttr->getName() == filter)
			{
				addAttribute(const_cast<Attribute*>(&*iAttr),
					iDefault->getBacktrackLevel(), attributes, attributeLookup);
			}
		}
	}
}

/**
 * Adds or updates the specified attribute in the lists.
 * @param attribute [in] the attribute to add (or update if it already exists).
 * @param backtrackLevel [in] the backtrack level to reach the attribute.
 * @param attributeList [out] the list of attributes where to add the new attribute.
 * @param attributeLookup [out] the map of existing attributes.
 */
void CConfig::addAttribute(Attribute* attribute,
	int backtrackLevel, AttributeList& attributeList,
	AttributeLookup& attributeLookup)
{
	// create a path node for the attribute to record its name and index.
	PathNode childNode(*attribute);
	// create a pseudo Result to store both the attribute pointer and the
	// backtrack level, and later associate them with a PathNode in the map.
	View::Result childSearchResult(attribute, backtrackLevel);

	// Look for the attribute in the map
	AttributeLookup::iterator iLookup = attributeLookup.find(childNode);
	// If the attribute does not already exists
	if(iLookup == attributeLookup.end())
	{
		// We add it to both lists
		attributeLookup.insert(make_pair(childNode, childSearchResult));
		attributeList.push_back(childNode);
	}
	else
	{
		// If it exists, we update the associated Result.
		// The Result references the nearest instance from the parent
		// object, which means: the instance that will override all
		// others.
		iLookup->second = childSearchResult;
	}
}


std::ostream& operator<<(std::ostream& stream, const CConfig& config)
{
	return stream << config.find().getValue();
}



///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig::Iterator class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor: creates an uninitialized iterator.
 */
CConfig::Iterator::Iterator()
	: filter(), parentConfig(), configs(), iConfig(), initialized(false)
{
}

/**
 * Constructor: initializes an iterator on the attributes of the specified object.
 * @param config the config on the object whose attributes will be enumerated.
 * @param filter enable to ignore attribute whose name does not match the filter.
 * @param includeDefault enable to include default index ('*') in iteration list.
 */
CConfig::Iterator::Iterator(const CConfig& config, const std::string& filter, bool includeDefault)
	: filter(filter), parentConfig(config), configs(), iConfig(), initialized(false), includeDefault(includeDefault)
{
}

/**
 * Copy constructor: creates a copy of an attribute's config iterator.
 * @param value the iterator to copy.
 */
CConfig::Iterator::Iterator(const Iterator& value)
	: filter(value.filter),
	parentConfig(value.parentConfig),
	configs(value.configs),
	iConfig(value.iConfig),
	initialized(value.initialized),
	includeDefault(value.includeDefault)
{
}

/**
 * Destructor.
 */
CConfig::Iterator::~Iterator()
{

}

/**
 * Assignment operator.
 * @param value the iterator to copy.
 */
CConfig::Iterator& CConfig::Iterator::operator =(const Iterator& value)
{
	filter = value.filter;
	parentConfig = value.parentConfig;
	configs = value.configs;
	iConfig = value.iConfig;
	initialized = value.initialized;
	includeDefault = value.includeDefault;
	return *this;
}

/**
 * Indicates whether there are more attributes left to enumerate.
 * @return Returns true if there are more attributes, false otherwise.
 * @throws CCfg::UnsupportedOperation if the underlying model attribute is not
 * an object.
 */
bool CConfig::Iterator::hasNext()
{
	if (not initialized)
	{
		initialize();
	}
	return iConfig + 1 != configs.end();
}

/**
 * Moves the iterator to the next attribute.
 * Any subsequent call to Iterator::get() will return the next attribute.
 * It is necessary to call Iterator::next() before calling Iterator::get() for
 * the first time.
 * hasNext must return true before calling this function, otherwise, it will
 * assert.
 * @return Returns the current attribute pointed by this iterator.
 * @throws CCfg::UnsupportedOperation if the underlying model attribute is not
 * an object.
 */
CConfig& CConfig::Iterator::next()
{
	if(hasNext())
	{
		++iConfig;
	}
	else
	{
		assert(false);
	}
	return *iConfig;
}

/**
 * Moves the iterator to the next attribute.
 * @see CCfg::CConfig::Iterator::next()
 */
CConfig::Iterator& CConfig::Iterator::operator ++()
{
	next();
	return *this;
}

/**
 * Returns the current attribute pointed by this iterator.
 * @see CCfg::CConfig::Iterator::get()
 */
CConfig& CConfig::Iterator::operator *()
{
	assert(iConfig !=  configs.end());
	return *iConfig;
}

/**
 * Allows to call functions on the current attribute pointed by this iterator.
 * @return Returns a pointer on the current attribute.
 */
CConfig* CConfig::Iterator::operator ->()
{
	assert(iConfig !=  configs.end());
	return &*iConfig;
}

/**
 * Returns the current attribute pointed by this iterator.
 * @return Returns the current attribute pointed by this iterator.
 */
CConfig& CConfig::Iterator::get()
{
	assert(iConfig !=  configs.end());
	return *iConfig;
}

/**
 * Returns the current attribute pointed by this iterator.
 * @see CCfg::CConfig::Iterator::get() const
 */
const CConfig& CConfig::Iterator::operator *() const
{
	assert(iConfig !=  configs.end());
	return *iConfig;
}

/**
 * Allows to call functions on the current attribute pointed by this iterator.
 * @return Returns a pointer on the current attribute.
 */
const CConfig* CConfig::Iterator::operator ->() const
{
	assert(iConfig !=  configs.end());
	return &*iConfig;
}

/**
 * Returns the current attribute pointed by this iterator.
 * @return Returns the current attribute pointed by this iterator.
 */
const CConfig& CConfig::Iterator::get() const
{
	assert(iConfig !=  configs.end());
	return *iConfig;
}

/**
 * Returns the number of element that will be enumerated.
 * @returns the number of element that will be enumerated.
 */
int CConfig::Iterator::getCount()
{
	if (not initialized)
	{
		initialize();
	}
	return configs.size();
}

void CConfig::Iterator::initialize()
{
	configs = parentConfig.getAttributes(filter, includeDefault);
	// The iteration will begin by a next(), therefore, we have to move one
	// step before the beginning.
	iConfig = configs.begin() - 1;
	initialized = true;
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig::NotFound class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param configPath the path of the config for which we couldn't find the
 * attribute.
 */
CConfig::NotFound::NotFound(const Path& configPath)
	: configPath(configPath)
{
}

/**
 * Copy constructor.
 * @param value the exception to copy.
 */
CConfig::NotFound::NotFound(const CConfig::NotFound& value)
	: Exception(value),	configPath(value.configPath)
{
}

/**
 * Destructor.
 */
CConfig::NotFound::~NotFound() throw()
{
}

/**
 * Assignment operator.
 * @param value the exception to copy
 */
CConfig::NotFound& CConfig::NotFound::operator =(
	const CConfig::NotFound& value)
{
	Exception::operator =(value);
	configPath = value.configPath;
	return *this;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string CConfig::NotFound::formatMessage() const
{
	ostringstream message;
	message << "The config element " << configPath << " cannot be resolved to a "
		"concrete attribute.";
	return message.str();
}

}
