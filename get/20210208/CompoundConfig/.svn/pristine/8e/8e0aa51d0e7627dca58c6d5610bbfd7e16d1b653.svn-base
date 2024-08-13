// =============================================================================
// Object.h
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

#ifndef CCfg_Object_h_INCLUDED
#define CCfg_Object_h_INCLUDED

#include "Attribute.h"
#include "Exception.h"
#include "Path.h"

#include <string>
#include <list>
#include <memory>

namespace CCfg
{

//forward declarations
class Visitor;

extern const char* DEFAULT_OBJ_INDEX;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a composite configuration object owning many children
 * sub-configurations.
 * These children are CCfg::Scalar and CCfg::Object and together with these
 * classes, they describe an object oriented, tree-like configuration.
 *
 * @author fchateau
 */
class Object : public Attribute
{
	std::string includeUrl;		///< the URL of the resource where the object's content comes from.
	std::list<Attribute*> attributes;	///< the children attributes of the object.

public:
	// construction/destruction
	Object(const std::string& name = "unnamed", const std::string& index = "");
	Object(const Object& value);
	~Object();
	std::auto_ptr<Object> clone(Recursivity recursivity = DEEP) const;

	// comparison
	bool equals(const Attribute& attr, MatchLevel level = IDENTICAL, Recursivity recursivity = DEEP) const;

	// visitation
	bool accept(Visitor& visitor) const;

	// iteration
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;

	// getters
	const std::string& getIncludeUrl() const;
	bool isDefaultConfig() const;

	// setters
	void setIncludeUrl(const std::string& value);
	void setDefaultConfig(bool value);

	// pseudo getters
	bool empty() const;
	int size() const;

	// insertion related methods
	void insert(iterator position, Attribute* value);
	void append(Attribute* attribute);
	void prepend(Attribute* attribute);
	using Attribute::append;

	// removal related methods
	std::auto_ptr<Attribute> remove(iterator attribute);
	void clear();
	using Attribute::remove;

	// high level config operations
	void merge(std::auto_ptr<Object> obj);
	void import(std::auto_ptr<Object> obj);

	// type information
	bool isObject() const;
	AttributeTypeId getType() const;

private:
	Object& operator =(const Object& value);		// assignment operator
	Attribute& operator =(const Attribute& value);	// overriden base class assignment
	Attribute& operator =(const String& value);
	Attribute& operator =(Int value);
	Attribute& operator =(Hex value);
	Attribute& operator =(Real value);
	Attribute& operator =(Bool value);
	Attribute& operator =(const char* value);
	Attribute& operator =(const Value& value);
	Attribute& operator =(short value);
	Attribute& operator =(long value);
	Attribute& operator =(unsigned short value);
	Attribute& operator =(unsigned long value);

	Object* cloneImpl(Recursivity recursivity = DEEP) const;
	void updateClassInfo();
};




///////////////////////////////////////////////////////////////////////////////
// Object inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the property isDefault.
 * @return Returns whether this configuration is a default configuration or not.
 */
inline bool Object::isDefaultConfig() const
{
	return getIndex() == DEFAULT_OBJ_INDEX;
}

/**
 * Setter of the property isDefault.
 * @param value the new value for the isDefault property.
 */
inline void Object::setDefaultConfig(bool value)
{
	if(value)
	{
		setIndex(DEFAULT_OBJ_INDEX);
	} else
	{
		setIndex("");
	}
}

}


#endif // CCfg_Object_h_INCLUDED
