// =============================================================================
// Attribute.h
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

#ifndef CCfg_Attribute_h_INCLUDED
#define CCfg_Attribute_h_INCLUDED

#include "Comment.h"
#include "Path.h"
#include "UnsupportedOperation.h"
#include "Value.h"
#include "TypeInfo.h"
#include <boost/shared_ptr.hpp>
#include <list>
#include <string>
#include <memory>
#include <iterator>

namespace CCfg
{

// forward declarations
class Object;
class Scalar;
class Visitor;
class ClassInfo;
class Range;
typedef std::auto_ptr<Range> RangePtr;

// enums
enum Recursivity
{
	SHALLOW,
	DEEP
};
enum MatchLevel
{
	EQUIVALENT, // requires structure, values, index, ranges, units to be the same.
	IDENTICAL // same as EQUIVALENT, but also requires comments, and include url to be the same.
};
enum InstanceType
{
	OBJECT_INST,
	SCALAR_INST
};

// typedefs
typedef boost::shared_ptr<ClassInfo> ClassInfoPtr;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

class Attribute;
typedef std::auto_ptr<Attribute> AttributePtr;

/**
 * This is the base class of the class hierarchy representing a configuration.
 * Configurations are made of CCfg::Objects and CCfg::Scalars. CCfg::Objects
 * own several Attributes that may be Objects or Scalars. Thus, these 3
 * classes implement a Composite pattern.
 *
 * @author fchateau
 */
class Attribute
{
protected:
	typedef std::list<Attribute*> AttrList;	///< the type of the list of attributes.

	Object* parentConfig;	///< the parent configuration object.
	std::string name;		///< the name of this configuration attribute.
	std::string index;		///< the index of this configuration attribute.
	Comment preComment;		///< the comment before this attribute.
	Comment postComment;	///< the comment after this attribute.
	ClassInfoPtr classInfo;	///< the shared metadata for members of the same "equivalence class".

public:
	// inner classes
	class iterator;
	class const_iterator;
	class NotFound;
	class NotAnAncestor;
	class IncompatibleClassInfo;

	// static methods
	static size_t computeAttributeCount(const Attribute& parent);

	// @{ Construction/destruction
	explicit Attribute(const std::string& name = "unnamed", const std::string& index = "");
	Attribute(const Attribute& value);
	virtual Attribute& operator =(const Attribute& value);
	virtual ~Attribute();
	// @}
	// cloning
	std::auto_ptr<Attribute> clone(Recursivity recursivity = DEEP) const;

	// high level operations
	void moveTo(Object& target);

	// @{ Navigation
	const Object* getIncludeRoot() const;
	Object* getRootConfig();
	// @}

	//@{ getters
	Object* getParent() const;
	const std::string& getName() const;
	const std::string& getIndex() const;
	Comment& getPreComment();
	Comment& getPostComment();
	const Comment& getPreComment() const;
	const Comment& getPostComment() const;
	const ClassInfoPtr& getClassInfo() const;
	//@}
	// @{ setters
	void setParent(Object *parent);
	void setName(const std::string& value);
	void setIndex(const std::string& value);
	void setPreComment(const Comment& value);
	void setPostComment(const Comment& value);
	void setClassInfo(const ClassInfoPtr& newClassInfo);
	// @}
	// @{ Path
	Path getPath(const Attribute* ancestor = 0) const;
	Path getClassPath() const;
	Path getShortClassPath() const;
	// @}
	// @{ Creation
	Attribute& append(InstanceType type, const std::string& name, const std::string& index = "");
	Attribute& require(InstanceType type, const std::string& name, const std::string& index = "");
	Attribute& require(InstanceType tailType, const Path& path);
	Attribute& requireScalar(const std::string& name, const std::string& index = "");
	Attribute& requireObject(const std::string& name, const std::string& index = "");
	Attribute& appendScalar(const std::string& name, const std::string& index = "");
	Attribute& appendObject(const std::string& name, const std::string& index = "");
	// @}
	// @{ Searches
	Attribute& operator ()(const char* name, const char* index = "");
	Attribute& operator ()(const char* name, int index);
	Attribute& operator ()(const std::string& name, const std::string& index = "");
	Attribute& operator ()(const std::string& name, int index);
	iterator find(const std::string& name, const std::string& index = "");
	iterator find(const std::string& name, int index);
	iterator search(const std::string& name, const std::string& index = "");
	iterator search(const std::string& name, int index);
	const_iterator find(const std::string& name, const std::string& index = "") const;
	const_iterator find(const std::string& name, int index) const;
	const_iterator search(const std::string& name, const std::string& index = "") const;
	const_iterator search(const std::string& name, int index) const;
	bool contains(const std::string& name, const std::string& index = "") const;
	// @}
	// @{ Iteration
	iterator getIterator();
	const_iterator getIterator() const;
	// @}
	// @{ Removal
	std::auto_ptr<Attribute> remove(Attribute* attr);
	std::auto_ptr<Attribute> remove(const std::string& name, const std::string& index = "");
	std::auto_ptr<Attribute> discard(const std::string& name, const std::string& index = "");
	// @}
	// @{ Conversion operators
	operator const String&() const;
	operator Int() const;
	operator Hex() const;
	operator Real() const;
	operator Bool() const;
	operator const char*() const;
	operator short() const;
	operator long() const;
	operator unsigned short() const;
	operator unsigned long() const;
	// @}
	// @{ Downcast methods
	Object& toObject();
	Scalar& toScalar();
	const Object& toObject() const;
	const Scalar& toScalar() const;
	// @}

	/*
	 * Pure interface.
	 */
	// @{ Misc
	virtual bool equals(const Attribute& attr, MatchLevel level = IDENTICAL, Recursivity recursivity = DEEP) const;
	virtual bool accept(Visitor& visitor) const = 0;
	// @}
	// @{ Type information
	virtual bool isObject() const = 0;
	virtual AttributeTypeId getType() const = 0;
	// @}

	/*
	 * Interface of CCfg::Object
	 */
	// @{ include
	virtual const std::string& getIncludeUrl() const;
	virtual void setIncludeUrl(const std::string& value);
	// @}
	// @{ Iteration
	virtual iterator begin();
	virtual iterator end();
	virtual const_iterator begin() const;
	virtual const_iterator end() const;
	// @}
	// @{ pseudo getters
	virtual int size() const;
	virtual bool empty() const;
	// @}
	// @{ Insertion
	virtual void insert(iterator position, Attribute* value);
	virtual void append(Attribute* subConfig);
	virtual void prepend(Attribute* subConfig);
	// @}
	// @{ Removal
	virtual std::auto_ptr<Attribute> remove(iterator attr);
	virtual void clear();
	// @}

	/*
	 * Interface of CCfg::Scalar
	 */
	// @{ Assignment operators
	virtual Attribute& operator =(const String& value) = 0;
	virtual Attribute& operator =(Int value) = 0;
	virtual Attribute& operator =(Hex value) = 0;
	virtual Attribute& operator =(Real value) = 0;
	virtual Attribute& operator =(Bool value) = 0;
	virtual Attribute& operator =(const char* value) = 0;
	virtual Attribute& operator =(const Value& value) = 0;
	virtual Attribute& operator =(short value) = 0;
	virtual Attribute& operator =(long value) = 0;
	virtual Attribute& operator =(unsigned short value) = 0;
	virtual Attribute& operator =(unsigned long value) = 0;
	// @}
	// @{ Unit/range
	virtual const std::string& getUnit() const;
	virtual Range* getRange() const;
	virtual void setUnit(const std::string& value);
	virtual void setRange(RangePtr range);
	virtual void setRange(const Range& range);
	// @}
	// @{ Getters
	virtual const Value&	getValue() const;
	virtual ValueTypeId		getValueTypeId() const;
	virtual const String&	getStringValue() const;
	virtual Int		getIntValue() const;
	virtual Hex		getHexValue() const;
	virtual Real	getRealValue() const;
	virtual Bool	getBoolValue() const;
	// @}
	// @{ Conversion getters
	virtual const String&	getValueAsString() const;
	virtual Int		getValueAsInt() const;
	virtual Hex		getValueAsHex() const;
	virtual Real	getValueAsReal() const;
	virtual Bool	getValueAsBool() const;
	// @}
	// @{ Setters
	virtual void setValue(const std::string& value);
	virtual void setValue(const Value& value);
	virtual void setValue(ValueTypeId type, const String& value);
	virtual void setStringValue(const String& value);

	virtual void setIntValue(Int value);
	virtual void setHexValue(Hex value);
	virtual void setRealValue(Real value);
	virtual void setBoolValue(Bool value);

	virtual void setIntValue(const String& value);
	virtual void setHexValue(const String& value);
	virtual void setRealValue(const String& value);
	virtual void setBoolValue(const String& value);
	// @}
protected:
	virtual void updateClassInfo() = 0;
	virtual Attribute* cloneImpl(Recursivity recursivity = DEEP) const = 0;

	UnsupportedOperation unsupported(const char* const opname) const;

	friend class Object;
};

std::ostream& operator<<(std::ostream& stream, const Attribute& attribute);


/**
 * The iterator on children attributes.
 */
class Attribute::iterator : public std::iterator<std::bidirectional_iterator_tag, Attribute>
{
	friend class Attribute::const_iterator;
	friend class Object;

	std::list<Attribute*>::iterator it;
public:
	// construction
	iterator() : it() {}
	iterator(const std::list<Attribute*>::iterator it) : it(it) {}
	iterator(const iterator& rval) : it(rval.it) {}
	~iterator() {}
	iterator& operator =(const iterator& rval) { it = rval.it; return *this; }

	// dereference
	Attribute& operator *() const { return **it; }
	Attribute* operator ->() const { return *it; }

	// incrementation
	iterator operator ++(int) { return iterator(it++); }
	iterator& operator ++() { ++it; return *this; }
	iterator operator --(int) { return iterator(it--); }
	iterator& operator --() { --it; return *this; }

	// comparison
	bool operator ==(const iterator& rval) { return it == rval.it; }
	bool operator !=(const iterator& rval) { return it != rval.it; }

	// extended methods (compared to standard iterators)
	// conversion
	operator Attribute*() { return *it; }
};



/**
 * The iterator on children attributes for read only purpose.
 */
class Attribute::const_iterator : public std::iterator<std::bidirectional_iterator_tag, const Attribute>
{
	std::list<Attribute*>::const_iterator it;
public:
	// construction
	const_iterator() : it() {}
	const_iterator(const std::list<Attribute*>::const_iterator it) : it(it) {}
	const_iterator(const std::list<Attribute*>::iterator it) : it(it) {}
	const_iterator(const Attribute::iterator& rval) : it(rval.it) {}

	const_iterator(const const_iterator& rval) : it(rval.it) {}
	~const_iterator() {}
	const_iterator& operator =(const const_iterator& rval) { it = rval.it; return *this; }

	// dereference
	const Attribute& operator *() const { return **it; }
	const Attribute* operator ->() const { return *it; }

	// incrementation
	const_iterator operator ++(int) { return const_iterator(it++); }
	const_iterator& operator ++() { ++it; return *this; }
	const_iterator operator --(int) { return const_iterator(it--); }
	const_iterator& operator --() { --it; return *this; }

	// comparison
	bool operator ==(const const_iterator& rval) { return it == rval.it; }
	bool operator !=(const const_iterator& rval) { return it != rval.it; }

	// extended methods (compared to standard iterators)
	// conversion
	operator const Attribute*() const { return *it; }
};

/**
 * This class describe the failure to find the requested attribute in a
 * configuration object.
 *
 * @author fchateau
 */
class Attribute::NotFound : public Exception
{
	Path objectPath;		///< the path of the object where the attribute was not found
	PathNode missingNode;	///< the name of the attribute that was supposed to exist

public:
	NotFound(const Attribute& attr, const std::string& name, const std::string& index);
	NotFound(const NotFound& value);
	~NotFound() throw();
	NotFound& operator =(const NotFound& value);

	const Path& getObjectPath() const { return objectPath; }
	const std::string& getName() const { return missingNode.name; }
	const std::string& getIndex() const { return missingNode.index; }

private:
	std::string formatMessage() const;
};

/**
 * This class reports a false assumption in the ascendency relation between
 * 2 attributes.
 * The error arise when an attribute, expected to be the ancestor of
 * another one in some function, has in fact no such relation.
 *
 * author fchateau
 */
class Attribute::NotAnAncestor : public Exception
{
	Path attributePath;	///< the path of the attribute for which we should have given a true ancestor.
	Path ancestorPath;	///< the path of the false ancestor attribute.

public:
	NotAnAncestor(const Attribute& attribute, const Attribute& ancestor);
	NotAnAncestor(const NotAnAncestor& value);
	~NotAnAncestor() throw();
	NotAnAncestor& operator =(const NotAnAncestor& value);

	const Path& getAttributePath() const;
	const Path& getAncestorPath() const;

private:
	std::string formatMessage() const;
};


/**
 * This class reports that an attribute joining a new equivalence class
 * already has a ClassInfo which is incompatible with the one of the class it joins.
 *
 * @author fchateau
 */
class Attribute::IncompatibleClassInfo : public Exception
{
	Path newComerPath;	///< the path of the new comer attribute in the equivalence class.
	ClassInfoPtr existingClassInfo;	///< the ClassInfo of equivalence class joined.
	ClassInfoPtr newComerClassInfo; ///< the ClassInfo of the former equivalence class of the new comer.

public:
	IncompatibleClassInfo(const Attribute* newComer, const Attribute* classMember);
	IncompatibleClassInfo(const IncompatibleClassInfo& rvalue);
	virtual ~IncompatibleClassInfo() throw();
	IncompatibleClassInfo& operator=(const IncompatibleClassInfo& right);

	const Path& getNewComerPath() const;
	const ClassInfo* getExistingClassInfo() const;
	const ClassInfo* getNewComerClassInfo() const;

private:
	std::string formatMessage() const;
};



///////////////////////////////////////////////////////////////////////////////
// CCfg::Attribute inline methods body
///////////////////////////////////////////////////////////////////////////////


/**
 * Getter of the "parentConfig" property.
 * @return Returns the parent Object config (null if this is a root)
 */
inline Object* Attribute::getParent() const
{
	return parentConfig;
}

/**
 * Setter of the "parentConfig" property.
 * @param value the new parent configuration.
 */
inline void Attribute::setParent(Object *value)
{
	parentConfig = value;
	updateClassInfo();
}

/**
 * Getter of the "name" property.
 * @return Returns the name of the configuration.
 */
inline const std::string& Attribute::getName() const
{
	return name;
}

/**
 * Getter of the "index" property.
 * @return Returns the index of the configuration.
 */
inline const std::string& Attribute::getIndex() const
{
	return index;
}

/**
 * Getter of the "preComment" property.
 * @return Returns the pre-comment for this configuration.
 */
inline Comment& Attribute::getPreComment()
{
	return preComment;
}

/**
 * Getter of the "postComment" property.
 * @return Returns the post-comment for this configuration.
 */
inline Comment& Attribute::getPostComment()
{
	return postComment;
}

/**
 * Getter of the "preComment" property.
 * @return Returns the pre-comment for this configuration.
 */
inline const Comment& Attribute::getPreComment() const
{
	return preComment;
}

/**
 * Getter of the "postComment" property.
 * @return Returns the post-comment for this configuration.
 */
inline const Comment& Attribute::getPostComment() const
{
	return postComment;
}

/**
 * Returns the shared metadata for members of the same "equivalence class".
 */
inline const ClassInfoPtr& Attribute::getClassInfo() const
{
	return classInfo;
}

/**
 * Setter of the "name" property.
 * @param value the new configuration name.
 */
inline void Attribute::setName(const std::string& value)
{
	name = value;
	updateClassInfo();
}

/**
 * Setter of the "index" property.
 * @param value the new configuration index.
 */
inline void Attribute::setIndex(const std::string& value)
{
	index = value;
}

/**
 * Setter of the "preComment" property.
 * @param value the new pre-comment.
 */
inline void Attribute::setPreComment(const Comment& value)
{
	preComment = value;
}

/**
 * Setter of the "postComment" property.
 * @param value the new post-comment.
 */
inline void Attribute::setPostComment(const Comment& value)
{
	postComment = value;
}

inline Attribute& Attribute::requireScalar(const std::string& name, const std::string& index)
{
	return require(SCALAR_INST, name, index);
}

inline Attribute& Attribute::requireObject(const std::string& name, const std::string& index)
{
	return require(OBJECT_INST, name, index);
}

inline Attribute& Attribute::appendScalar(const std::string& name, const std::string& index)
{
	return append(SCALAR_INST, name, index);
}

inline Attribute& Attribute::appendObject(const std::string& name, const std::string& index)
{
	return append(OBJECT_INST, name, index);
}

/**
 * Converts this attribute to its scalar string value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a string.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator const String&() const
{
	return getStringValue();
}

/**
 * Converts this attribute to its scalar integer value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an integer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator Int() const
{
	return getIntValue();
}

/**
 * Converts this attribute to its scalar hexadecimal value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an hexadecimal number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator Hex() const
{
	return getHexValue();
}

/**
 * Converts this attribute to its scalar real value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a real number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator Real() const
{
	return getRealValue();
}

/**
 * Converts this attribute to its scalar boolean value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a boolean.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator Bool() const
{
	return getBoolValue();
}

/**
 * Converts this attribute to its scalar string value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a string.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator const char*() const
{
	return getStringValue().c_str();
}

/**
 * Converts this attribute to its scalar integer value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an integer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator short() const
{
	return getIntValue();
}

/**
 * Converts this attribute to its scalar integer value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an integer.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator long() const
{
	return getIntValue();
}

/**
 * Converts this attribute to its scalar hexadecimal value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an hexadecimal number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator unsigned short() const
{
	return getHexValue();
}

/**
 * Converts this attribute to its scalar hexadecimal value, if possible.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an hexadecimal number.
 * @throws UnsupportedOperation if this CCfg::Attribute is not a CCfg::Scalar.
 */
inline Attribute::operator unsigned long() const
{
	return getHexValue();
}

/**
 * Stream insertion operator.
 */
inline std::ostream& operator<<(std::ostream& stream, const Attribute& attribute)
{
	return stream << attribute.getValue();
}

}

#endif // CCfg_Attribute_h_INCLUDED
