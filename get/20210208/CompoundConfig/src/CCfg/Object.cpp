// =============================================================================
// Object.cpp
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

#include "Object.h"
#include "Scalar.h"
#include "Visitor.h"
#include "Attribute.h"
#include <sstream>
#include <algorithm>
#include <functional>

#include <assert.h>

using std::string;
using std::ostringstream;
using std::unary_function;
using std::count_if;

namespace CCfg
{
//static members initialization
const char* DEFAULT_OBJ_INDEX = "*";


///////////////////////////////////////////////////////////////////////////////
// Object class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: creates an object config with the specified name and index.
 * @param name the name of the new config object.
 * @param index the index of the new config object.
 */
Object::Object(const std::string& name, const std::string& index)
	: Attribute(name, index)
{
}

/**
 * Copy constructor: do a shallow copy of an object (do not copy children).
 * @param value the config object to copy.
 */
Object::Object(const Object& value)
	: Attribute(value), includeUrl(value.includeUrl)
{
}

/**
 * Destructor: recursively invokes delete on any children attribute.
 */
Object::~Object()
{
	clear();
}

/**
 * Makes a clone of this config object.
 * @note
 * The clone is not included as a child of the parent of its model.
 * @param recursivity allows to select between deep copy or shallow copy.
 * @return Returns a newly allocated copy of this config object.
 */
Object* Object::cloneImpl(Recursivity recursivity) const
{
	Object *c = new Object(*this);
	if(recursivity == DEEP)
	{
		for(const_iterator iAttribute = begin(), endAttribute = end();
			iAttribute != endAttribute; ++iAttribute)
		{
			c->append(iAttribute->clone(recursivity).release());
		}
	}
	return c;
}

/*
 *
 */
void Object::updateClassInfo()
{
	for(iterator iAttribute = begin(), endAttribute = end();
		iAttribute != endAttribute; ++iAttribute)
	{
		iAttribute->updateClassInfo();
	}
}

/**
 * Makes a clone of this config object.
 * @note
 * The clone is not included as a child of the parent of its model.
 * @param recursivity allows to select between deep copy or shallow copy.
 * @return Returns a newly allocated copy of this object.
 */
std::auto_ptr<Object> Object::clone(Recursivity recursivity) const
{
	return std::auto_ptr<Object>(cloneImpl(recursivity));
}

/*
 * (see inherited doc)
 */
bool Object::equals(const Attribute& attr, MatchLevel level, Recursivity recursivity) const
{
	const Object *obj = dynamic_cast<const Object*>(&attr);
	if((obj != 0 and	// check if attr is also an Object.
		Attribute::equals(attr, level, recursivity) and	// check if bases classes are equal
		level == EQUIVALENT) or
		(level == IDENTICAL and getIncludeUrl() == obj->getIncludeUrl()))
	{
		if(recursivity == DEEP)
		{
			const_iterator iLeft = begin(), iRight = obj->begin();
			while(iLeft != end() and iRight != obj->end())
			{
				// the first inequality make the comparison fails
				if(not iLeft->equals(*iRight, level, recursivity))
				{
					return false;
				}
				 ++iLeft, ++iRight;
			}
			// every attributes are equal.

			// if one object's content has not been fully enumerated
			// one is longer that the other.
			if(iLeft != end() or iRight != obj->end())
			{
				return false;
			}
			else
			{
				return true;
			}
		}
		else
		{
			return true;
		}
	}
	return false;
}

/**
 * Implement visitor support.
 * Iterates recursively on every children configuration attributes and call
 * the visitor method corresponding to its type.
 * @param visitor the visitor object.
 * @return Returns false if the visit was interrupted by the visitor (it happens
 * when one of its visit method returns false).
 */
bool Object::accept(Visitor& visitor) const
{
	if(visitor.onObjectEnter(*this))
	{
		// iterates on each attributes
		for(const_iterator iAttribute = begin(), endAttribute = end();
			iAttribute != endAttribute; ++iAttribute)
		{
			if(not iAttribute->accept(visitor))
			{
				// false indicates that visitation should ignore next siblings.
				break;
			}
		}

	}
	return visitor.onObjectExit(*this);
}

/**
 * Gives the starting point for an iteration loop on children.
 * @return returns an iterator on the first children of this object.
 */
Object::iterator Object::begin()
{
	return iterator(attributes.begin());
}

/**
 * Gives the end point for an iteration loop on children.
 * @return Returns an iterator pointing over the last children of this attribute.
 */
Object::iterator Object::end()
{
	return iterator(attributes.end());
}

/**
 * Const-correct version of Object::begin().
 */
Object::const_iterator Object::begin() const
{
	return const_iterator(attributes.begin());
}

/**
 * Const-correct version of Object::end().
 */
Object::const_iterator Object::end() const
{
	return const_iterator(attributes.end());
}

/**
 * Getter of the "includeUrl" property.
 * @return Returns the URL of the included resource if there is one (and
 * an empty string if there is none).
 */
const std::string& Object::getIncludeUrl() const
{
	return includeUrl;
}

/**
 * Setter of the "includeUrl" property.
 * @param value the URL of the new included resource.
 */
void Object::setIncludeUrl(const std::string& value)
{
	includeUrl = value;
}

/**
 * Indicates whether this attribute has children or not.
 * @return Returns true if this attribute has no children.
 */
bool Object::empty() const
{
	return attributes.empty();
}

/**
 * Gives the number of children config attributes.
 * This function is not constant time but O(n), so use it sparingly.
 * In particular in situations where "empty" might be sufficient.
 * @return Returns the number of children config attributes.
 */
int Object::size() const
{
	return attributes.size();
}

/**
 * Inserts a config attribute as a children of this object, at a specific location.
 * @param position an iterator on the attribute before which it will be inserted.
 * @param subConfig the config attribute to insert.
 * @throws UnsupportedOperation if this Attribute is not a CCfg::Object.
 */
void Object::insert(iterator position, Attribute* subConfig)
{
	assert(subConfig != 0);
	assert(subConfig->getParent() == 0);

	attributes.insert(position.it, subConfig);
	subConfig->setParent(this);
}

/**
 * Inserts a config attribute as a children of this object at the end of the list.
 * @param attribute the config attribute to insert.
 */
void Object::append(Attribute* attribute)
{
	assert(attribute != 0);
	assert(attribute->getParent() == 0);

	attributes.push_back(attribute);
	attribute->setParent(this);
}

/**
 * Inserts a config attribute as a children of this object at the beginning of the list.
 * @param attribute the config attribute to insert.
 */
void Object::prepend(Attribute* attribute)
{
	assert(attribute != 0);
	assert(attribute->getParent() == 0);

	attributes.push_front(attribute);
	attribute->setParent(this);
}

/**
 * Removes the specified attribute and its sub-tree from its parent.
 * @param attr an iterator on the root attribute of the sub-tree to remove.
 * @return Returns a pointer to the root of the removed sub-tree.
 * The pointer ownership is given back to the caller through an auto pointer.
 */
std::auto_ptr<Attribute> Object::remove(iterator attr)
{
	assert(attr != attributes.end());

	std::auto_ptr<Attribute> removed(attr);
	attributes.erase(attr.it);
	removed->setParent(0);
	return removed;
}

/**
 * Deletes all children of this attribute.
 * @throws UnsupportedOperation if this Attribute is not a CCfg::Object.
 */
void Object::clear()
{
	for(iterator iAttr = attributes.begin(), endAttr = attributes.end();
		iAttr != endAttr; ++iAttr)
	{
		delete iAttr;
	}
	attributes.clear();
}

/**
 * Merges the content of an object with this one.
 * @note
 * The content of the source object will be moved to the target. The source
 * object must be a standalone configuration or must have been previously
 * removed because the method takes ownership of the pointer and destroy it
 * before returning.
 * @param obj the object whose content is to be merged.
 */
void Object::merge(std::auto_ptr<Object> obj)
{
	iterator itNext = obj->begin();
	iterator it;
	while(itNext != obj->end())
	{
		it = itNext++;
		it->moveTo(*this);
		// iterate on itNext because moveTo deletes items in the list while
		// iterating, and it would be unsafe to delete the item pointed
		// by the current iterator.
	}

}

/**
 * Import the content of an object in this one.
 * @note
 * The content of the source object will be moved to the target. The source
 * object must be a standalone configuration or must have been previously
 * removed because the method takes ownership of the pointer and destroy it
 * before returning.
 * The content of the destination object will be erased beforehand.
 */
void Object::import(std::auto_ptr<Object> obj)
{
	clear();
	merge(obj);
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(const Attribute& /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(const String& /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(Int /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(Hex /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(Real /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(Bool /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(const char* /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(const Value& /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(short /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(long /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(unsigned short /*value*/)
{
	throw unsupported("operator =");
}

/**
 * Disabled assignment operator. Reserved for scalars.
 * @throws UnsupportedOperation
 */
Attribute& Object::operator =(unsigned long /*value*/)
{
	throw unsupported("operator =");
}

/*
 * (see inherited doc)
 */
bool Object::isObject() const
{
	return true;
}

/*
 * (see inherited doc)
 */
AttributeTypeId Object::getType() const
{
	return OBJECT;
}




}

