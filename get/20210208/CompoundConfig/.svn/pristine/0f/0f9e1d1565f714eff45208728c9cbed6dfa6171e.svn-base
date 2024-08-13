// =============================================================================
// CConfig.h
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

#ifndef CCfg_CConfig_h_INCLUDED
#define CCfg_CConfig_h_INCLUDED

#include "View/Result.h"
#include "Object.h"
#include "Path.h"
#include "Exception.h"
#include "TypeInfo.h"
#include <vector>
#include <map>

namespace CCfg
{

// forward declarations
class Value;
class Comment;
class CConfig;

namespace View
{
class InheritedIterator;
class InheritersIterator;
class ClassIterator;
}

// enums
/// A list of configs generated when enumerating the content of the parent config.
typedef std::vector<CConfig> ConfigList;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class manipulates the configuration model using the "default value"
 * semantics.
 *
 * @author fchateau
 */
class CConfig
{
	/**
	 * An attribute of the configuration to manipulate.
	 * It is the attribute from where the path resolution starts.
	 */
	Attribute* baseAttribute;
	/**
	 * The current path of the virtual attribute represented by the config.
	 * It is altered each time the user calls operator() or find.
	 * This is a transient value.
	 */
	mutable Path path;
	/**
	 * The index of the path node in "path" corresponding to "baseAttribute".
	 * It is needed to reset the path after resolution.
	 */
	int baseDepth;
	/**
	 * The number of backtrackings when searching the base attribute from the
	 * root, using the path.
	 */
	int backtrackLevel;

public:
	class Iterator;
	class NotFound;

	// @{ construction/destruction
	CConfig();
	CConfig(Attribute* baseAttribute);
	CConfig(Attribute* baseAttribute, const Path& path, int backtrackLevel);
	CConfig(const CConfig& config);
	~CConfig();
	CConfig& operator =(const CConfig& value);
	// @}

	// @{ getters/setters
	Attribute* getBaseAttribute() const;
	void setBaseAttribute(Attribute* newBase);
	int getBacktrackLevel() const;
	// @}

	// @{ Resolution
	Attribute& create(InstanceType lastInstType);
	Attribute* search() const;
	Attribute& find() const;
	View::Result resolve() const;
	bool exists() const;
	// @}

	 // @{ path management
	const Path& getPath() const;
	bool hasBasePath() const;
	CConfig& operator ()(const char* name, const char* index = "");
	CConfig& operator ()(const char* name, int index);
	CConfig& operator ()(const std::string& name, const std::string& index = "");
	CConfig& operator ()(const std::string& name, int index);
	CConfig& at(const char* relativePath);
	CConfig& append(const std::string& name, const std::string& index = "");
	CConfig& append(const std::string& name, int index);
	CConfig& appendPath(const Path& relativePath);
	CConfig& appendPath(const std::string& relativePath);
	// @}

	// @{ iteration
	Iterator iterate(std::string nameFilter = "", bool includeDefault = true) const;
	ConfigList getAttributes(const std::string& nameFilter = "", bool includeDefault = true) const;
	View::InheritedIterator		getInheritedAttributes() const;
	View::InheritersIterator	getInheriterAttributes() const;
	View::ClassIterator			getClassInstances() const;
	// @}

private:
	void resetPath() const;

public:
	/*
	 * Delegation interface to CCfg::Attribute
	 */
	// @{ getters
	const std::string& getName() const;
	const std::string& getIndex() const;
	Comment& getPreComment();
	Comment& getPostComment();
	const Comment& getPreComment() const;
	const Comment& getPostComment() const;
	const std::string& getIncludeUrl() const;
	// @}

	// @{ setters
	void setName(const std::string& value);
	void setIndex(const std::string& value);
	void setPreComment(const Comment& value);
	void setPostComment(const Comment& value);
	void setIncludeUrl(const std::string& value);
	// @}

	// @{ scalar assignment operators
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
	template<typename T> Attribute& assign(T value);
	// @}

	// @{ scalar conversion operators
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

	/*
	 * Delegation interface to CCfg::Scalar
	 */
	// @{ Unit/range
	const std::string& getUnit() const;
	Range* getRange() const;
	void setUnit(const std::string& value);
	void setRange(RangePtr range);
	void setRange(const Range& range);
	// @}
	// @{ getters
	const Value&		getValue() const;
	ValueTypeId			getValueTypeId() const;

	const String& getStringValue () const;
	Int           getIntValue    () const;
	Hex           getHexValue    () const;
	Real          getRealValue   () const;
	Bool          getBoolValue   () const;

	const String& getStringValue (const String& defaultValue) const;
	Int           getIntValue    (Int defaultValue) const;
	Hex           getHexValue    (Hex defaultValue) const;
	Real          getRealValue   (Real defaultValue) const;
	Bool          getBoolValue   (Bool defaultValue) const;
	// @}

	// @{ conversion getters
	const String&	getValueAsString() const;
	Int		getValueAsInt() const;
	Hex		getValueAsHex() const;
	Real	getValueAsReal() const;
	Bool	getValueAsBool() const;
	// @}

	// @{ setters
	void setValue(const std::string& value);
	void setValue(const Value& value);
	void setStringValue(const String& value);

	void setIntValue(Int value);
	void setHexValue(Hex value);
	void setRealValue(Real value);
	void setBoolValue(Bool value);

	void setIntValue(const String& value);
	void setHexValue(const String& value);
	void setRealValue(const String& value);
	void setBoolValue(const String& value);
	// @}

private:
	/// the list of config resolution results of each inherited object of the
	/// object whose content is to be enumerated.
	typedef std::vector<View::Result> InheritedList;

	/// the list of attributes of a config. Only name and index is needed, hence,
	/// we use PathNode objects.
	typedef std::vector<PathNode> AttributeList;

	/// a lookup table to know whether an attribute has already been found.
	/// It is associated with a Result object because we will need the attribute
	/// pointer and the backtrack level to build a config for this attribute.
	typedef std::map<PathNode, View::Result> AttributeLookup;


	void listInherited(InheritedList& defaults) const;

	static void listAttributes(const InheritedList& defaults,
		const std::string& filter,
		AttributeList& attributes,
		AttributeLookup& attributeLookup);

	static void addAttribute(Attribute* attribute,
		int backtrackLevel,
		AttributeList& attributes,
		AttributeLookup& attributeLookup);
};

std::ostream& operator<<(std::ostream& stream, const CConfig& config);

/**
 * This class allows to enumerate the whole list of attribute of an object,
 * even those inherited from a default object.
 *
 * @author fchateau
 */
class CConfig::Iterator
{
	std::string filter;
	CConfig parentConfig;	///< the config on the object whose attributes will be enumerated.
	ConfigList configs;		///< the list of configs created for each attribute.
	ConfigList::iterator iConfig;	///< the iterator on the current config.
	bool initialized; ///< whether the iterator has been initialized
	bool includeDefault; ///< Whether default index is included in iteration

public:
	// construction/destruction
	Iterator();
	Iterator(const CConfig& config, const std::string& filter = "", bool includeDefault = true);
	Iterator(const Iterator& value);
	virtual ~Iterator();
	Iterator& operator =(const Iterator& value);

	// iteration
	bool hasNext();
	CConfig& next();
	Iterator& operator ++();

	// dereference
	CConfig& operator *();
	CConfig* operator ->();
	CConfig& get();
	const CConfig& operator *() const;
	const CConfig* operator ->() const;
	const CConfig& get() const;

	// other
	int getCount();

private:
	void initialize();
};



/**
 * This class describe the failure to find the requested attribute in a
 * configuration.
 *
 * @author fchateau
 */
class CConfig::NotFound : public Exception
{
	Path configPath;	///< the path of the attribute that cannot be found

public:
	NotFound(const Path& configPath);
	NotFound(const NotFound& value);
	~NotFound() throw();
	NotFound& operator =(const NotFound& value);

	const Path& getConfigPath() const;

private:
	std::string formatMessage() const;
};




///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the property "baseAttribute".
 * @return Returns the attribute of the configuration where the config starts.
 */
inline Attribute* CConfig::getBaseAttribute() const
{
	return baseAttribute;
}

/**
 * Getter of the property "path".
 * @return Returns the path of the virtual attribute represented by the config.
 */
inline const Path& CConfig::getPath() const
{
	return path;
}

/**
 * Returns whether the path is the base path or if it has been appended.
 */
inline bool CConfig::hasBasePath() const
{
	// if the base depth is equal to the path size - 1, it means nothing has been appended to the base path
	return path.size()-1 == baseDepth;
}

/**
 * Getter of the property "backtrackLevel".
 * @return Returns the number of backtrackings when searching the base
 * attribute from the root, using the path.
 */
inline int CConfig::getBacktrackLevel() const
{
	return backtrackLevel;
}

/**
 * Shortcut for the append method.
 * @see append(const std::string&, const std::string&)
 */
inline CConfig& CConfig::operator ()(const std::string& name, const std::string& index)
{
	return append(name, index);
}

/**
 * Shortcut for the append method.
 * @see append(const std::string&, int)
 */
inline CConfig& CConfig::operator ()(const std::string& name, int index)
{
	return append(name, index);
}

/**
 * Returns whether the attribute referenced by the current path exists.
 */
inline bool CConfig::exists() const
{
	return search() != 0;
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig::NotFound inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the property "configPath".
 * @return Returns the path that the config was not able to translate into a
 * concrete attribute.
 */
inline const Path& CConfig::NotFound::getConfigPath() const
{
	return configPath;
}

}

#endif // CCfg_CConfig_h_INCLUDED
