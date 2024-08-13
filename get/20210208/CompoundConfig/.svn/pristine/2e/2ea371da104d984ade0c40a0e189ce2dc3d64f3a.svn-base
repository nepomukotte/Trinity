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

#include "View/ClassIterator.h"
#include "Attribute.h"
#include "Object.h"
#include "Scalar.h"
#include "Comment.h"
#include "Visitor.h"
#include "ClassInfo.h"
#include "Range.h"
#include "Utils/Io/PointerValueStreamer.hpp"
#include <sstream>
#include <cassert>

using std::string;
using std::ostringstream;
using Utils::Io::writeValue;

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Attribute class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: creates an attribute with the specified name and index.
 * @param name the name of the new config attribute.
 * @param index the index of the new config attribute.
 */
Attribute::Attribute(const std::string& name, const std::string& index)
	: parentConfig(NULL), name(name), index(index),
	preComment(), postComment(), classInfo()
{
}

/**
 * Copy constructor: do a shallow copy of an attribute.
 * It does not copy children, nor integrate the new attribute in the source
 * tree.
 * @param value the config attribute to copy.
 */
Attribute::Attribute(const Attribute& value)
	: parentConfig(NULL), name(value.name), index(value.index),
	preComment(value.preComment), postComment(value.postComment),
	classInfo()
{
	if(value.classInfo.get() != 0)
	{
		std::auto_ptr<ClassInfo> rightClass = value.classInfo->clone();
		classInfo = rightClass;
	}
}

/**
 * Destructor
 */
Attribute::~Attribute()
{
}

/**
 * Assignment operator: do a shallow copy of a scalar attribute.
 * @note
 * This operator is virtual to allow assignment between scalars manipulated
 * from their base classes.
 * @param value the config attribute to copy.
 * @throws UnsupportedOperation if the left and right values are not of the same type.
 */
Attribute& Attribute::operator =(const Attribute& value)
{
	// the parent is not changed
	setName(value.getName());
	setIndex(value.getIndex());
	setPreComment(value.getPreComment());
	setPostComment(value.getPostComment());
	// the class info should not be copied
	return *this;
}

/**
 * Makes a clone of this config attribute.
 * @note
 * The clone is not included as a child of the parent of its model.
 * @param recursivity allows to select between deep copy or shallow copy.
 * @return Returns a newly allocated copy of this attribute.
 */
std::auto_ptr<Attribute> Attribute::clone(Recursivity recursivity) const
{
	return std::auto_ptr<Attribute>(cloneImpl(recursivity));
}

/**
 * Compares this attribute to another for equality.
 * @param attr the comparison attribute.
 * @param level indicates the level of similarity required.
 * @param recursivity indicates whether comparison should be deep or shallow.
 * A shallow comparison ignores children attributes, and a deep comparison
 * requires both subtrees to be the same.
 */
bool Attribute::equals(const Attribute& attr, MatchLevel level, Recursivity /*recursivity*/) const
{
	return getName() == attr.getName() and
		getIndex() == attr.getIndex() and
		(level == EQUIVALENT or
		(level == IDENTICAL and
		getPreComment() == attr.getPreComment() and
		getPostComment() == attr.getPostComment()));
}


/**
 * Moves an attribute (and its sub-tree) to another object.
 * @param target the object where it will be added as a child.
 */
void Attribute::moveTo(Object& target)
{
	try
	{
		Object* parent = getParent();
		if(parent != 0)
		{
			parent->remove(this).release();
			// do no keep ownership of the pointer because it is the this
			// pointer and deleting it here would be an undefined operation.
		}
		target.append(this);
	}
	catch(NotAnAncestor& e)
	{
		// this should never happen. otherwise, it's a bug.
		assert(false);
	}
}

/**
 * Returns the first parent object having an include Url.
 */
const Object* Attribute::getIncludeRoot() const
{
	const Attribute* current = this;
	if(current->getParent() != 0)
	{
		do
		{
			current = current->getParent();
		}
		while(current->getParent() != 0 and current->getIncludeUrl().empty());
	}
	else
	{
		assert(current->isObject());
	}
	return static_cast<const Object*>(current);
}

/**
 * Returns the root of the current attribute.
 * @return Returns the root of the current attribute.
 */
Object* Attribute::getRootConfig()
{
	Attribute* currentAttr = this;
	while(currentAttr->getParent() != 0)
	{
		currentAttr = currentAttr->getParent();
	}
	return static_cast<Object*>(currentAttr);
}

/**
 * Defines the ClassInfo of every member of the equivalence class of this attribute.
 */
void Attribute::setClassInfo(const ClassInfoPtr& newClassInfo)
{
	View::ClassIterator it(*this);
	while(it.hasNext())
	{
		it.next()->classInfo = newClassInfo;
	}
}

/**
 * Returns the qualified path of the attribute relative to an ancestor.
 * @param ancestor the ancestor attribute from which the path is expressed.
 * If null, the full path is computed (including the root attribute).
 * @returns Returns the qualified path of the attribute in its configuration.
 * @throws Attribute::NotAncestor if the specified ancestor attribute is not an
 * ancestor of the attribute for which we compute the path.
 */
Path Attribute::getPath(const Attribute* ancestor) const
{
	Path path;
	const Attribute *current = this;

	while(current != 0 and current != ancestor)
	{
		path.prepend(*current);
		current = current->getParent();
	}

	// if we found the ancestor at the end of the search.
	if(current == ancestor)
	{
		return path;
	} else
	{
		throw NotAnAncestor(*this, *ancestor);
		// dereferencing is ok: if ancestor == 0,
		// "current" always equals 0 at the end
	}
}

/**
 * Returns the path of the attribute's class.
 */
Path Attribute::getClassPath() const
{
	Path path;
	const Attribute *current = this;

	while(current != 0)
	{
		path.prepend(PathNode(current->getName(), "*"));
		current = current->getParent();
	}
	return path;
}

/**
 * Returns the path of the attribute's class in its short form (without wildcards).
 */
Path Attribute::getShortClassPath() const
{
	Path path;
	const Attribute *current = this;

	while(current != 0)
	{
		path.prepend(PathNode(current->getName(), string()));
		current = current->getParent();
	}
	return path;
}

/**
 * Creates a new attribute and appends it to the list of children attribute.
 */
Attribute& Attribute::append(InstanceType type, const std::string& name, const std::string& index)
{
	Attribute* attr = (type == OBJECT_INST) ?
			static_cast<Attribute*>(new Object(name, index)) :
			static_cast<Attribute*>(new Scalar(name, index));
	append(attr);
	return *attr;
}

/**
 * Requires the existence of a children attribute, creating it if necessary.
 */
Attribute& Attribute::require(InstanceType type, const std::string& name, const std::string& index)
{
	iterator iAttr = search(name, index);
	if(iAttr == end())
	{
		return append(type, name, index);
	}
	else
	{
		return *iAttr;
	}
}

/**
 * Requires the existence of a chain of children attributes, creating some of
 * them if necessary.
 */
Attribute& Attribute::require(InstanceType tailType, const Path& path)
{
	Attribute* currentAttribute = this;
	Attribute* nextAttribute;
	for(Path::const_iterator iNode = path.begin(), endNode = path.end();
		iNode != endNode; ++iNode)
	{
		iterator iAttr = currentAttribute->search(iNode->name, iNode->index);
		if(iAttr == currentAttribute->end())
		{
			if(iNode+1 == endNode)
			{
				nextAttribute = &currentAttribute->require(tailType, iNode->name, iNode->index);
			}
			else
			{
				nextAttribute = &currentAttribute->requireObject(iNode->name, iNode->index);
			}
		}
		else
		{
			nextAttribute = iAttr;
		}
		currentAttribute = nextAttribute;
	}
	return *currentAttribute;
}

/**
 * Returns the attribute having the specified name and index.
 */
Attribute& Attribute::operator ()(const char* name, const char* index)
{
	return *find(name, index);
}

/**
 * Returns the attribute having the specified name and index.
 */
Attribute& Attribute::operator ()(const char* name, int index)
{
	return *find(name, index);
}

/**
 * Returns the attribute having the specified name and index.
 */
Attribute& Attribute::operator ()(const std::string& name, const std::string& index)
{
	return *find(name, index);
}

/**
 * Returns the attribute having the specified name and index.
 */
Attribute& Attribute::operator ()(const std::string& name, int index)
{
	return *find(name, index);
}

/**
 * Finds the children attribute having the specified name and index.
 * \note
 * Find and search distinguish themselves by the fact that find is supposed
 * to always success. So it must throw an exception otherwise. On the other
 * hand, when failure is part of the normal flow, one should use "search".
 * @param name the name of the attribute to find.
 * @param index the index of the attribute to find.
 * @return Returns an iterator on the matching children attribute.
 * @throws Attribute::NotFound if there is no children attribute matching the
 * criteria.
 */
Attribute::iterator Attribute::find(const std::string& name, const std::string& index)
{
	iterator result = search(name, index);
	if(result == end())
	{
		throw NotFound(*this, name, index);
	} else
	{
		return result;
	}
}

/**
 * Const-correct version of Object::find(const std::string&, const std::string&).
 */
Attribute::const_iterator Attribute::find(const std::string& name, const std::string& index) const
{
	const_iterator result = search(name, index);
	if(result == end())
	{
		throw NotFound(*this, name, index);
	} else
	{
		return result;
	}
}

/**
 * Finds the children attribute having the specified name and index.
 * \note
 * Find and search distinguish themselves by the fact that find is supposed
 * to always success. So it must throw an exception otherwise. On the other
 * hand, when failure is part of the normal flow, one should use "search".
 * @param name the name of the attribute to find.
 * @param index the index of the attribute to find.
 * @return Returns an iterator on the matching children attribute.
 * @throws Attribute::NotFound if there is no children attribute matching the
 * criteria.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
Attribute::iterator Attribute::find(const std::string& name, int index)
{
	iterator result = search(name, index);
	if(result == end())
	{
		ostringstream indexStr;
		indexStr << index;
		throw NotFound(*this, name, indexStr.str());
	} else
	{
		return result;
	}
}

/**
 * Const-correct version of Object::find(const std::string&, int).
 */
Attribute::const_iterator Attribute::find(const std::string& name, int index) const
{
	const_iterator result = search(name, index);
	if(result == end())
	{
		ostringstream indexStr;
		indexStr << index;
		throw NotFound(*this, name, indexStr.str());
	} else
	{
		return result;
	}
}

/**
 * Look for a children attribute having the specified name and index.
 * \note
 * Find and search distinguish themselves by the fact that find is supposed
 * to always success. So it must throw an exception otherwise. On the other
 * hand, when failure is part of the normal flow, one should use "search".
 * @param name the name of the attribute to search.
 * @param index the index of the attribute to search.
 * @return Returns an iterator to the matching children attribute.
 * If it doesn't exist, it is equal to Object::end().
 */
Attribute::iterator Attribute::search(const std::string& name, const std::string& index)
{
	iterator iAttribute = begin(), endAttribute = end();
	for(; iAttribute != endAttribute; ++iAttribute)
	{
		if(iAttribute->getName() == name and iAttribute->getIndex() == index)
		{
			return iAttribute;
		}
	}
	return endAttribute;
}

/**
 * Const-correct version of Object::search(const std::string&, const std::string&).
 */
Attribute::const_iterator Attribute::search(const std::string& name, const std::string& index) const
{
	const_iterator iAttribute = begin(), endAttribute = end();
	for(; iAttribute != endAttribute; ++iAttribute)
	{
		if(iAttribute->getName() == name and iAttribute->getIndex() == index)
		{
			return iAttribute;
		}
	}
	return endAttribute;
}

/**
 * Look for a children attribute having the specified name and index.
 * \note
 * Find and search distinguish themselves by the fact that find is supposed
 * to always success. So it must throw an exception otherwise. On the other
 * hand, when failure is part of the normal flow, one should use "search".
 * @param name the name of the attribute to search.
 * @param index the index of the attribute to search.
 * @return Returns an iterator to the matching children attribute.
 * If it doesn't exist, it is equal to Object::end().
 */
Attribute::iterator Attribute::search(const std::string& name, int index)
{
	ostringstream indexStr;
	indexStr << index;
	return search(name, indexStr.str());
}

/**
 * Const-correct version of Object::search(const std::string&, const std::string&).
 */
Attribute::const_iterator Attribute::search(const std::string& name, int index) const
{
	ostringstream indexStr;
	indexStr << index;
	return search(name, indexStr.str());
}

/**
 * Returns whether this attribute contains an attribute of the specified name and index.
 */
bool Attribute::contains(const std::string& name, const std::string& index) const
{
	return search(name, index) != end();
}

/**
 * Returns the iterator on this attribute in its parent.
 */
Attribute::iterator Attribute::getIterator()
{
	assert(parentConfig != 0);
	iterator it = parentConfig->begin();
	iterator end = parentConfig->end();
	while(it != end and it != this)
	{
		++it;
	}
	assert(it == this);
	return it;
}

/**
 * Returns the iterator on this attribute in its parent.
 */
Attribute::const_iterator Attribute::getIterator() const
{
	assert(parentConfig != 0);
	const_iterator it = parentConfig->begin();
	const_iterator end = parentConfig->end();
	while(it != end and it != this)
	{
		++it;
	}
	assert(it == this);
	return it;
}


/**
 * Removes the specified attribute and its sub-tree from its parent.
 * @param attr the root of the sub-tree to remove.
 * @return Returns a pointer to the root of the removed sub-tree.
 * The pointer ownership is given back to the caller through an auto pointer.
 * @throws Object::NotAnAncestor reports this object is not the ancestor of the
 * attribute.
 */
std::auto_ptr<Attribute> Attribute::remove(Attribute* attr)
{
	assert(attr != 0);
	if(attr->getParent() != this)
	{
		throw NotAnAncestor(*this, *attr);
	}
	iterator it = attr->getIterator();
	return remove(it);
}

/**
 * Removes the matching attribute and its sub-tree from its parent.
 * @param name the name of the attribute to remove.
 * @param index the index of the attribute to remove.
 * @return Returns the root of the removed sub-tree.
 * @throws Attribute::NotFound if there is no children attribute matching the criteria.
 */
std::auto_ptr<Attribute> Attribute::remove(const std::string& name, const std::string& index)
{
	return remove(find(name, index));
}

/**
 * Removes
 */
std::auto_ptr<Attribute> Attribute::discard(const std::string& name, const std::string& index)
{
	iterator iAttr = search(name, index);
	if(iAttr != end())
	{
		return remove(iAttr);
	}
	else
	{
		return std::auto_ptr<Attribute>();
	}
}

/**
 * Downcasts a config attribute to a config object.
 * This is a more readable alternative to downcast with dynamic_cast.
 * @return Returns a reference on the config object.
 * @throws std::bad_cast if this CCfg::Attribute is not a CCfg::Object.
 */
Object& Attribute::toObject()
{
	return dynamic_cast<Object&>(*this);
}

/**
 * Downcasts a config attribute to a config scalar.
 * This is a more readable alternative to downcast with dynamic_cast.
 * @return Returns a reference on the config scalar.
 * @throws std::bad_cast if this attribute is not a scalar.
 */
Scalar& Attribute::toScalar()
{
	return dynamic_cast<Scalar&>(*this);
}

/**
 * Downcasts a config attribute to a config object.
 * This is a more readable alternative to downcast with dynamic_cast.
 * @return Returns a reference on the config object.
 * @throws std::bad_cast if this CCfg::Attribute is not a CCfg::Object.
 */
const Object& Attribute::toObject() const
{
	return dynamic_cast<const Object&>(*this);
}

/**
 * Downcasts a config attribute to a config scalar.
 * This is a more readable alternative to downcast with dynamic_cast.
 * @return Returns a reference on the config scalar.
 * @throws std::bad_cast if this attribute is not a scalar.
 */
const Scalar& Attribute::toScalar() const
{
	return dynamic_cast<const Scalar&>(*this);
}

/**
 * @fn bool Attribute::isObject() const
 * Returns whether this attribute is an object or not.
 */

/**
 * @fn AttributeTypeId Attribute::getType() const
 * Returns the type of this attribute.
 */

/**
 * Recursively computes the number of attributes in a configuration subtree.
 * @param parent the root of the subtree.
 */
size_t Attribute::computeAttributeCount(const Attribute& parent)
{
	struct CounterVisitor : Visitor
	{
		size_t count;
		CounterVisitor() : Visitor(), count(0) {}
		~CounterVisitor() {}
		bool onObjectEnter(const Object&) { ++count; return true; }
		bool onObjectExit(const Object&) { return true; }
		bool onScalar(const Scalar&) { ++count; return true; }
	} v;
	parent.accept(v);
	return v.count;
}






UnsupportedOperation Attribute::unsupported(const char* const opname) const
{
	return UnsupportedOperation(getPath().toString(), typeid(*this), opname);
}

/**
 * Getter of the "includeUrl" property.
 * @return Returns the URL of the included resource if there is one (and
 * an empty string if there is none).
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
const std::string& Attribute::getIncludeUrl() const
{
	throw unsupported("getIncludeUrl");
}

/**
 * Setter of the "includeUrl" property.
 * @param value the URL of the new included resource.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
void Attribute::setIncludeUrl(const std::string& /*value*/)
{
	throw unsupported("setIncludeUrl");
}

/**
 * Gives the starting point for an iteration loop on children.
 * @return returns an iterator on the first children of this attribute.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
Attribute::iterator Attribute::begin()
{
	throw unsupported("begin");
}

/**
 * Gives the end point for an iteration loop on children.
 * @return Returns an iterator pointing over the last children of this attribute.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
Attribute::iterator Attribute::end()
{
	throw unsupported("end");
}

/**
 * Const-correct version of Attribute::begin().
 */
Attribute::const_iterator Attribute::begin() const
{
	throw unsupported("begin");
}

/**
 * Const-correct version of Attribute::end().
 */
Attribute::const_iterator Attribute::end() const
{
	throw unsupported("end");
}

/**
 * Indicates whether this attribute has children or not.
 * @return Returns true if this attribute has no children.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
bool Attribute::empty() const
{
	throw unsupported("empty");
}

/**
 * Gives the number of children config attributes.
 * This function is not constant time but O(n), so use it sparingly,
 * particularily in situations where "empty" might be sufficient.
 * @return Returns the number of children config attributes.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
int Attribute::size() const
{
	throw unsupported("size");
}

/**
 * Inserts a config attribute as a children of this one, at a specific location.
 * @param position an iterator on the attribute before which it will be inserted.
 * @param attribute the config attribute to insert.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
void Attribute::insert(Attribute::iterator /*position*/, Attribute* /*attribute*/)
{
	throw unsupported("insert");
}

/**
 * Inserts a config attribute as a children of this one at the end of the list.
 * @param attribute the config attribute to insert.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
void Attribute::append(Attribute* /*attribute*/)
{
	throw unsupported("append");
}

/**
 * Inserts a config attribute as a children of this one at the beginning of the list.
 * @param attribute the config attribute to insert.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
void Attribute::prepend(Attribute* /*attribute*/)
{
	throw unsupported("prepend");
}


/**
 * Removes the specified attribute and its sub-tree from its parent.
 * @param attr an iterator on the root attribute of the sub-tree to remove.
 * @return Returns a pointer to the root of the removed sub-tree.
 * The pointer ownership is given back to the caller through an auto pointer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
std::auto_ptr<Attribute> Attribute::remove(iterator /*attr*/)
{
	throw unsupported("remove");
}

/**
 * Deletes all children of this attribute.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Object.
 */
void Attribute::clear()
{
	throw unsupported("clear");
}



/**
 * Getter of the property "unit".
 * @return Returns the physical unit of the value.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
const std::string& Attribute::getUnit() const
{
	throw unsupported("getUnit");
}

/**
 * Returns the range of a Scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Range* Attribute::getRange() const
{
	throw unsupported("getRange");
}

/**
 * Setter of the propery "unit".
 * @param value the new value for this property.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setUnit(const std::string& /*value*/)
{
	throw unsupported("setUnit");
}

/**
 * Defines the range of the equivalence class of this Scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setRange(RangePtr /*range*/)
{
	throw unsupported("setRange");
}

/**
 * Defines the range of the equivalence class of this Scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setRange(const Range& /*range*/)
{
	throw unsupported("setRange");
}

/**
 * Getter of the property "val".
 * @return Returns the scalar value.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
const Value& Attribute::getValue() const
{
	throw unsupported("getValue");
}

/**
 * Getter of the property "valueType".
 * @return Returns the type of the value in the scalar (integer, real, etc...)
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
ValueTypeId Attribute::getValueTypeId() const
{
	throw unsupported("getValueTypeId");
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a string.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
const String& Attribute::getStringValue() const
{
	throw unsupported("getStringValue");
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an integer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Int Attribute::getIntValue() const
{
	throw unsupported("getIntValue");
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an hexadecimal number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Hex Attribute::getHexValue() const
{
	throw unsupported("getHexValue");
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a real number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Real Attribute::getRealValue() const
{
	throw unsupported("getRealValue");
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a boolean.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Bool Attribute::getBoolValue() const
{
	throw unsupported("getBoolValue");
}

/**
 * Getter the scalar value, with conversion is necessary.
 * @return Returns the scalar value converted to string.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
const String& Attribute::getValueAsString() const
{
	throw unsupported("getValueAsString");
}

/**
 * Getter the scalar value, with conversion is necessary.
 * @return Returns the scalar value converted to integer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Int	Attribute::getValueAsInt() const
{
	throw unsupported("getValueAsInt");
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to hexadecimal.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Hex	Attribute::getValueAsHex() const
{
	throw unsupported("getValueAsHex");
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to real.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Real Attribute::getValueAsReal() const
{
	throw unsupported("getValueAsReal");
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to boolean.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
Bool Attribute::getValueAsBool() const
{
	throw unsupported("getValueAsBool");
}

/**
 * Defines the scalar value.
 * @param value the new value for this property.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setValue(const std::string& /*value*/)
{
	throw unsupported("setValue");
}

/**
 * Defines the scalar value.
 * @param value the new value for this property.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setValue(const Value& /*value*/)
{
	throw unsupported("setValue");
}

/**
 * Defines the scalar value.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setValue(ValueTypeId /*type*/, const String& /*value*/)
{
	throw unsupported("setValue");
}

/**
 * Defines the scalar value and set its type to string.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setStringValue(const String& /*value*/)
{
	throw unsupported("setStringValue");
}

/**
 * Defines the scalar value and set its type to integer.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setIntValue(Int /*value*/)
{
	throw unsupported("setIntValue");
}

/**
 * Defines the scalar value and set its type to hexadecimal.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setHexValue(Hex /*value*/)
{
	throw unsupported("setHexValue");
}

/**
 * Defines the scalar value and set its type to integer.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setRealValue(Real /*value*/)
{
	throw unsupported("setRealValue");
}

/**
 * Defines the scalar value and set its type to boolean.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setBoolValue(Bool /*value*/)
{
	throw unsupported("setBoolValue");
}

/**
 * Defines the scalar value and set its type to integer.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setIntValue(const String& /*value*/)
{
	throw unsupported("setIntValue");
}

/**
 * Defines the scalar value and set its type to hexadecimal.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setHexValue(const String& /*value*/)
{
	throw unsupported("setHexValue");
}

/**
 * Defines the scalar value and set its type to integer.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setRealValue(const String& /*value*/)
{
	throw unsupported("setRealValue");
}

/**
 * Defines the scalar value and set its type to boolean.
 * @param value the new value for the scalar.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
void Attribute::setBoolValue(const String& /*value*/)
{
	throw unsupported("setBoolValue");
}

/**
 * @fn Attribute* Attribute::cloneImpl(Recursivity recursivity) const
 * Makes a clone of this config attribute.
 * @param recursivity allows to select between deep copy or shallow copy.
 * @return Returns a newly allocated copy of this config object.
 */

/**
 * @fn void Attribute::updateClassInfo()
 * Recursively updates links to equivalence classes, starting from *this.
 * This function is automatically called when the name or the parent object of
 * an attribute changes. It ensures the consistency of equivalence classes and
 * always keep them up to date.
 * @throws Attribute::IncompatibleClassInfo
 */






///////////////////////////////////////////////////////////////////////////////
// Attribute::NotFound class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param attr the object where the attribute was not found
 * @param name the name of the attribute that was supposed to exist
 * @param index the index of the attribute that was supposed to exist
 */
Attribute::NotFound::NotFound(const Attribute& attr,
	const std::string& name, const std::string& index)
	: objectPath(attr.getPath()), missingNode(name, index)
{
}

/**
 * Copy constructor.
 * @param value the exception to copy.
 */
Attribute::NotFound::NotFound(
	const Attribute::NotFound& value)
	: Exception(value),	objectPath(value.objectPath),
	missingNode(value.missingNode)
{
}

/**
 * Destructor.
 */
Attribute::NotFound::~NotFound() throw()
{
}

/**
 * Assignment operator.
 * @param value the exception to copy
 */
Object::Attribute::NotFound& Attribute::NotFound::operator =(
	const Attribute::NotFound& value)
{
	Exception::operator =(value);
	objectPath = value.objectPath;
	missingNode = value.missingNode;
	return *this;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Attribute::NotFound::formatMessage() const
{
	ostringstream message;
	message << "Cannot find the attribute " << missingNode
		<< " in the object: " << objectPath;
	return message.str();
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Attribute::NotAnAncestor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param attribute the attribute for which we should have given a true ancestor.
 * @param ancestor the false ancestor attribute.
 */
Attribute::NotAnAncestor::NotAnAncestor(const Attribute& attribute, const Attribute& ancestor)
	: Exception(),
	attributePath(attribute.getPath()),
	ancestorPath(ancestor.getPath())
{
}

/**
 * Copy constructor.
 * @param value the exception to copy.
 */
Attribute::NotAnAncestor::NotAnAncestor(const NotAnAncestor& value)
	: Exception(value),
	attributePath(value.getAttributePath()),
	ancestorPath(value.getAncestorPath())
{
}

/**
 * Destructor.
 */
Attribute::NotAnAncestor::~NotAnAncestor() throw()
{
}

/**
 * Assignment operator.
 * @param value the exception to copy
 */
Attribute::NotAnAncestor& Attribute::NotAnAncestor::operator =(const NotAnAncestor& value)
{
	Exception::operator =(value);
	attributePath = value.getAttributePath();
	ancestorPath = value.getAncestorPath();
	return *this;
}

/**
 * Getter of the "attributePath" property.
 * @return Returns the path of the attribute for which we should have given a
 * true ancestor.
 */
const Path& Attribute::NotAnAncestor::getAttributePath() const
{
	return attributePath;
}

/**
 * Getter of the "ancestorPath" property.
 * @return Returns the path of the false ancestor attribute.
 */
const Path& Attribute::NotAnAncestor::getAncestorPath() const
{
	return ancestorPath;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Attribute::NotAnAncestor::formatMessage() const
{
	ostringstream message;
	message << "Error in getPath method: the attribute " << attributePath
		<< " is not an ancestor of " << ancestorPath;
	return message.str();
}



////////////////////////////////////////////////////////////////////////////////
// CCfg::Attribute::IncompatibleClassInfo class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param newComer the attribute that has just joined the incompatible equivalence class.
 * @param classMember a member of the joined equivalence class.
 */
Attribute::IncompatibleClassInfo::IncompatibleClassInfo(const Attribute* newComer,
	const Attribute* classMember)
	: Exception(),
	newComerPath(newComer->getPath()),
	existingClassInfo(classMember->getClassInfo()),
	newComerClassInfo(newComer->getClassInfo())
{
}

/**
 * Copy constructor.
 */
Attribute::IncompatibleClassInfo::IncompatibleClassInfo(const IncompatibleClassInfo& rvalue)
	: Exception(rvalue),
	newComerPath(rvalue.newComerPath),
	existingClassInfo(rvalue.existingClassInfo),
	newComerClassInfo(rvalue.newComerClassInfo)
{
}

/**
 * Destructor.
 */
Attribute::IncompatibleClassInfo::~IncompatibleClassInfo() throw()
{
}

/**
 * Assignment operator.
 */
Attribute::IncompatibleClassInfo&
Attribute::IncompatibleClassInfo::operator=(const IncompatibleClassInfo& rvalue)
{
	Exception::operator=(rvalue);
	newComerPath = rvalue.newComerPath;
	existingClassInfo = rvalue.existingClassInfo;
	newComerClassInfo = rvalue.newComerClassInfo;
	return *this;
}

/**
 * Returns the path of the new comer attribute in the equivalence class.
 */
const Path& Attribute::IncompatibleClassInfo::getNewComerPath() const
{
	return newComerPath;
}

/**
 * Returns the ClassInfo of equivalence class joined.
 */
const ClassInfo* Attribute::IncompatibleClassInfo::getExistingClassInfo() const
{
	return existingClassInfo.get();
}

/**
 * Returns the ClassInfo of the former equivalence class of the new comer.
 */
const ClassInfo* Attribute::IncompatibleClassInfo::getNewComerClassInfo() const
{
	return newComerClassInfo.get();
}

/*
 *
 */
std::string Attribute::IncompatibleClassInfo::formatMessage() const
{
	ostringstream message;
	message << "The informations about the configuration class of the attribute "
		<< newComerPath << " which are: " << writeValue(newComerClassInfo.get())
		<< " are incompatible with those of the class it has joined: "
		<< writeValue(existingClassInfo.get());
	return message.str();
}

}
