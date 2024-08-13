// =============================================================================
// BinaryFormat.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Kirell Benzi
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Io_BinaryFormat_h_INCLUDED
#define CCfg_Io_BinaryFormat_h_INCLUDED

#include "../AttributeClass.h"
#include "CCfg/Io/Format.h"
#include "CCfg/Attribute.h"
#include "CCfg/TypeInfo.h"
#include <boost/cstdint.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <map>
#include <vector>
#include <set>
#include <iostream>

// Make sure host endianness has been properly defined (typically in config.h file).
// If not, compilation is aborted.
#if not defined(WORDS_BIG_ENDIAN) and not defined(WORDS_LITTLE_ENDIAN)
#error "Endianness has not been defined!"
#endif

namespace CCfg
{

// Forward declarations
class Attribute;
class Object;
class Scalar;

//Typedefs
typedef uint64_t HexWord;
typedef uint64_t RealWord;
typedef uint32_t TableSize;
typedef uint32_t IntWord;
typedef uint32_t StringId;
typedef uint32_t AttributeClassId;
typedef uint16_t ChildCount;
typedef uint16_t VersionNumber;
typedef uint8_t  BoolWord;
typedef uint8_t  Byte;
typedef uint8_t  AttributeTypeWord;// enum: Object, Int, Real, Hex,Bool, String
typedef uint8_t  RangeCodeWord; //NO_RANGE, INTERVAL, ENUM
typedef uint8_t  RangeTypeWord; //IncludeBound, ExcludeBound, Unbounded

typedef std::vector<Byte> ByteBuffer;
typedef std::map<std::string, StringId> StringToIdMap;
typedef std::map<std::string, AttributeClassId> AttributeClassToIdMap;
typedef std::map<StringId,std::string> IdToStringMap;
typedef boost::ptr_map<AttributeClassId, AttributeClass> IdToAttributeClassMap;

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////
/**
 * This class provide methods to map the configuration semantics in ByteCode.
 *
 * @author kbenzi
 * */

class BinaryFormat : public Format
{
public:
	/// Inner classes
	// @{
	class SemanticError;
	// @}

	BinaryFormat();

	///Virtual function implementation
	std::auto_ptr<CCfg::Object> parse(Storage& storage, IncludeHandler& handler);

	///Virtual function implementation
	//@{
	void format(const Object *config, Storage& storage, IncludeHandler& handler);
	void format(const Attribute* attr);
	virtual std::string getConfigName(Storage& storage);
	//@}

	//To manipulate the Version number
	VersionNumber getVersion() const { return version; }
	void setVersion(VersionNumber);

	//Get the buffer
	ByteBuffer getIoBuffer() { return ioBuffer; }

	Object* parse();
	//Virtual method implementation
	std::auto_ptr<CCfg::Object> parse(Storage& storage);

//***************** FILE READER methods ***********//
	void readBinaryFromFile(std::string url);

private:
//***************** FORMAT methods ****************//
	static void serialize2cstring(ByteBuffer& buff, const std::string& str);

	//Recursive call
	void recursiveFormat(const Attribute* attr);

	///@name Serialize Range
	//@{
	void serializeRange(ByteBuffer& buff, const ValueRange* vRange);
	template<typename T, typename FormatT>
	void serializeBoundaries(const ValueRange* vRange);
	template<typename T, typename FormatT>
	void serializeEnum(const ValueRange* vRange);
	//@}

	///@name Serializer
	//@{
	template<typename TargetT, typename SourceT>
	void serializePrimitive(ByteBuffer& buff, SourceT value);
	void serializeAttributeClass(const AttributeClass* attr);
	void serializeAttribute(const Attribute* attr, AttributeClassId attrClassId);
	//@}

	StringId processString(const std::string& str);
	void processAttribute(const Attribute* attr);
	AttributeClassId processAttributeClass(const AttributeClass* attrClass);

	void repack(); //Gather all the buffers
	void cleanUp(); //clean all the buffers

	///@name Range Parsers
	//@{
	template<typename TargetT, typename SourceT>
	ValueRangeInterval<TargetT>* parseRangeInterval(ByteBuffer& buff, size_t& entryPoint);
	template<typename TargetT, typename SourceT>
	ValueRangeEnum<TargetT>* parseRangeEnum(ByteBuffer& buff, size_t& entryPoint);
	template<typename WordT>
	RangePtr makeInterval(ValueRange* valueRange);
	template<typename WordT>
	RangePtr makeEnum(ValueRange* valueRange);
	//@}

	///@name Parsers
	//@{
	template<typename TargetT, typename SourceT>
	TargetT parsePrimitive(ByteBuffer& buff, size_t& entryPoint);
	std::auto_ptr<AttributeClass> parseAttributeClass(size_t& currentPos);
	void parseStringMap(TableSize stringTNum,size_t& currentPos );
	void parseAttributeClassMap(TableSize attrClassNum,size_t& currentPos);
	Attribute* parseAttribute(ByteBuffer& buff,size_t& currentPos);
	//@}

	//Reformat a ValueRange into a Range format
	RangePtr reaRange(ValueRange* v);


private:

	//***** Maps for formatting *****//
	/// Map used to store all the strings in the configuration and to get an Id
	StringToIdMap stringToIdMap;
	/// Map used to store the name of each configClass to get an Id
	AttributeClassToIdMap attributeClassToIdMap;
	//***** Maps for parsing *****//
	/// Map used to get a String from an Id
	IdToStringMap idToStringMap;
	/// Map used to get AttributeClass from an Id
	IdToAttributeClassMap idToAttributeClassMap;

	//Byte Buffers
	ByteBuffer stringTableBuffer;
	ByteBuffer attributeClassBuffer;
	ByteBuffer attributeBuffer;
	ByteBuffer ioBuffer; // BYTE CODE BUFFER
	VersionNumber version;
	TableSize stringCount;
	TableSize attributeClassCount;

};

/**
 * Reports that the binary stream does not respect the semantics of our configuration.
 */
class BinaryFormat::SemanticError : public Format::SemanticError
{
public:
	SemanticError();
	~SemanticError() throw();
private:
	std::string formatMessage() const;
};

template<>
inline void BinaryFormat::serializePrimitive<StringId,String>(ByteBuffer& buff, String value)
{
	StringId strId = processString(value);
	serializePrimitive<StringId,StringId>(buff,strId);
}

/**
 * This template serializes a SourceT value into a ByteBuffer after transformation.
 * The real value written is a TargetT type which is a kind of cast on the SourceT value.
 * Depending on the host machine endianness the serialization is different.
 * @param buff ByteBuffer where the value is serialized
 * @param value SourceT value to serialize
 * @tparam SourceT original type of the value
 * @tparam TargetT final value after transformation
 */
template<typename TargetT, typename SourceT>
void BinaryFormat::serializePrimitive(ByteBuffer& buff, SourceT value)
{
	const int sizeT = static_cast<int>(sizeof(TargetT));
	const int sizeS = static_cast<int>(sizeof(SourceT));

	int deltaBytes = sizeS - sizeT;
	Byte* dataPtr = reinterpret_cast<Byte*>(&value);

#ifdef WORDS_BIG_ENDIAN // if host is big endian
	while (deltaBytes < 0)
	{
		buff.push_back(static_cast<Byte>(0)); // Zero-pad because target is wider than source
		deltaBytes++;
	}
	while (deltaBytes < sizeS)
	{
		buff.push_back(dataPtr[deltaBytes]);
		deltaBytes++;
	}
#else
	while (deltaBytes < 0)
	{
		buff.push_back(0); // Zero-pad because target is wider than source
		deltaBytes++;
	}
	while (deltaBytes < sizeS)
	{
		buff.push_back(dataPtr[sizeS - 1 - deltaBytes]);
		deltaBytes++;
	}
#endif
}


template<>
inline void BinaryFormat::serializeBoundaries<String,StringId>(const ValueRange* vRange)
{
	const ValueRangeInterval<String>& vRInterval = *(static_cast<const ValueRangeInterval<String>*> (vRange));
	for(int i=0; i<2; i++)
	{
		serializePrimitive<RangeTypeWord, Interval::BoundMode>(attributeClassBuffer, vRInterval[i].boundMode);
		StringId strId = processString(vRInterval[i].value);
		serializePrimitive<StringId, StringId>(attributeClassBuffer,strId);
	}
}

/**
 * Template function used to serialize an ValueRangeInterval of an AttributeClass.
 * First it serializes a type of Interval (included bounds, etc.).
 * Then the T-Interval type
 * @param vRange Range in the AttributeClass to serialize
 * @tparam T Type of the Range
 * @tparam FormatT Type the range need to be casted to
 */
template<typename T, typename FormatT>
void BinaryFormat::serializeBoundaries(const ValueRange* vRange)
{
	const ValueRangeInterval<T>& vRInterval = *(static_cast<const ValueRangeInterval<T>*> (vRange));
	for(int i=0; i<2; i++)
	{
		serializePrimitive<RangeTypeWord, Interval::BoundMode>(attributeClassBuffer, vRInterval[i].boundMode);
		serializePrimitive<FormatT, T>(attributeClassBuffer,vRInterval[i].value);
	}
}

template<>
inline void BinaryFormat::serializeEnum<String,StringId>(const ValueRange* vRange)
{
	const ValueRangeEnum<String>& vREnum = *(static_cast<const ValueRangeEnum<String>*> (vRange));
	std::set<String>::const_iterator it;

	serializePrimitive<TableSize,TableSize>(attributeClassBuffer,vREnum.values.size());

	for( it = vREnum.values.begin(); it != vREnum.values.end() ; ++it)
	{
		StringId strId = processString(*it);
		serializePrimitive<StringId, StringId>(attributeClassBuffer,strId);
	}
}

/**
 * Template function used to serialize an ValueRangeEnum of an AttributeClass.
 * @param vRange Range in the AttributeClass to serialize
 * @tparam T Type of the Range
 * @tparam FormatT Type the range need to be casted to
 */
template<typename T, typename FormatT>
void BinaryFormat::serializeEnum(const ValueRange* vRange)
{
	const ValueRangeEnum<T>& vREnum = *(static_cast<const ValueRangeEnum<T>*> (vRange));
	typename std::set<T>::const_iterator it;

	//Serialize the number of enum
	serializePrimitive<TableSize,TableSize>(attributeClassBuffer,vREnum.values.size());

	for( it = vREnum.values.begin(); it != vREnum.values.end() ; ++it)
	{
		serializePrimitive<FormatT, T>(attributeClassBuffer,*it);
	}
}

template<>
inline String BinaryFormat::parsePrimitive<String,StringId>(ByteBuffer& buff, size_t& entryPoint)
{
	StringId strId = parsePrimitive<StringId,StringId>(buff,entryPoint);
	IdToStringMap::iterator itStr = idToStringMap.find(strId); //finds the right iterator on string
	if(itStr == idToStringMap.end())
	{
		throw SemanticError();
	}
	return itStr->second; // return the string
}

/**
 * Gets a value from a Byte buffer and returns a variable of TargetT type.
 * The method manages the endianness.
 * @param buff ByteBuffer where the value is parsed from
 * @param entryPoint entrypoint in the buffer
 * @tparam SourceT original type of the value in the buffer
 * @tparam TargetT final type of the value returned
 * @return Value from the buffer after transformation
 */
template<typename TargetT, typename SourceT>
TargetT BinaryFormat::parsePrimitive(ByteBuffer& buff, size_t& entryPoint)
{
	const int sizeT = static_cast<int>(sizeof(TargetT));
	const int sizeS = static_cast<int>(sizeof(SourceT));

	int deltaBytes =  sizeT - sizeS;
	TargetT  targetVal;
	Byte*    targetPtr = reinterpret_cast<Byte*>(&targetVal);

#ifdef WORDS_BIG_ENDIAN // if host is big endian

	// Skipping buffer MSB if buffer word is wider than host word
	while (deltaBytes < 0)
	{
		entryPoint++;
		deltaBytes++;
	}

	// Padding host word if wider than buffer word
	for (int i=0; i < deltaBytes; ++i)
	{
		targetPtr[i] = 0;
	}

	// Copying buffer word LSB into host word LSB
	while (deltaBytes < sizeT)
	{
		targetPtr[deltaBytes] = buff.at(entryPoint++);
		deltaBytes++;
	}

#else

	// Skipping buffer MSB if buffer word is wider than host word
	while (deltaBytes < 0)
	{
		entryPoint++;
		deltaBytes++;
	}

	// Padding host word if wider than buffer word
	for (int i = 0; i < deltaBytes; ++i)
	{
		targetPtr[sizeT-1-i] = static_cast<TargetT>(0);
	}

	// Copying buffer word LSB into host word LSB
	for (int i = sizeT-deltaBytes-1 ; i >= 0 ; i--)
	{
		targetPtr[i] = buff.at(entryPoint++);
	}

#endif

	return targetVal;
}


template<>
inline ValueRangeInterval<String>* BinaryFormat::parseRangeInterval<String,StringId>(ByteBuffer& buff, size_t& entryPoint)
{
	ValueRangeInterval<String>* rangeInter = new ValueRangeInterval<String>();
	for(int i = 0; i < 2; ++i)
	{
	rangeInter->bounds[i].boundMode = parsePrimitive<Interval::BoundMode,RangeTypeWord>(buff,entryPoint);
	rangeInter->bounds[i].value = parsePrimitive<String,StringId>(buff,entryPoint);
	}
	return rangeInter;
}

/**
 * parseRangeInterval gets from the byteBuffer an interval and returns it.
 * First the boundMode and the values
 * @param buff ByteBuffer where the interval is parsed from
 * @param entryPoint entrypoint in the buffer
 * @tparam SourceT original type of the values in the buffer
 * @tparam TargetT final type of the values returned
 * @return ValueRangeInterval of TargetT type
 */
template<typename TargetT, typename SourceT>
ValueRangeInterval<TargetT>* BinaryFormat::parseRangeInterval(ByteBuffer& buff, size_t& entryPoint)
{
	ValueRangeInterval<TargetT>* rangeInter = new ValueRangeInterval<TargetT>();
	for(int i = 0; i < 2; ++i)
	{
		rangeInter->bounds[i].boundMode = parsePrimitive<Interval::BoundMode,RangeTypeWord>(buff,entryPoint);
		rangeInter->bounds[i].value = parsePrimitive<TargetT,SourceT>(buff,entryPoint);
	}
	return rangeInter;
}


template<>
inline ValueRangeEnum<String>* BinaryFormat::parseRangeEnum<String,StringId>(ByteBuffer& buff, size_t& entryPoint)
{
	ValueRangeEnum<String>* rangeEnum = new ValueRangeEnum<String>();
	TableSize num = parsePrimitive<TableSize,TableSize>(buff,entryPoint); //gets the numbers of enums

	for (TableSize i = 0; i < num; ++i) //for each enum
		rangeEnum->values.insert(parsePrimitive<String,StringId>(buff,entryPoint));

	return rangeEnum;
}

/**
 * Gets the number of enums and fills the enum attribute
 * obtained from the buffer.
 * @param buff ByteBuffer where the interval is parsed from
 * @param entryPoint entrypoint in the buffer
 * @tparam SourceT original type of the values parsed from the buffer
 * @tparam TargetT final type of the values inserted in the object
 * @return a pointer on a filled ValueRangeEnum<T>
 */
template<typename TargetT, typename SourceT>
ValueRangeEnum<TargetT>* BinaryFormat::parseRangeEnum(ByteBuffer& buff, size_t& entryPoint)
{
	ValueRangeEnum<TargetT>* rangeEnum = new ValueRangeEnum<TargetT>();
	TableSize num = parsePrimitive<TableSize,TableSize>(buff,entryPoint); //gets the numbers of enums

	for (TableSize i = 0; i< num; ++i) //for each enum
	{
		rangeEnum->values.insert(parsePrimitive<TargetT,SourceT>(buff,entryPoint));
	}
	return rangeEnum;
}

/**
 * Creates an Interval object from a ValueRange object
 * @param valueRange Original valueRange
 * @tparam WordT Interval Type
 * @return Range pointer on the new interval
 */
template<typename WordT>
RangePtr BinaryFormat::makeInterval(ValueRange* valueRange)
{
	ValueRangeInterval<WordT>* vRI;
	vRI = dynamic_cast< ValueRangeInterval<WordT>* >(valueRange);
	return RangePtr
	(
		new Interval
		(
			vRI->bounds[0].value,
			vRI->bounds[1].value,
			vRI->bounds[0].boundMode,
			vRI->bounds[1].boundMode
		)
	);
}

/**
 * Creates an Enum object from a ValueRange object
 * @param valueRange Original valueRange
 * @tparam WordT Enum Type
 * @return Range pointer on the new enum
 */
template<typename WordT>
RangePtr BinaryFormat::makeEnum(ValueRange* valueRange)
{
	ValueRangeEnum<WordT>* vRE;
	vRE = dynamic_cast< ValueRangeEnum<WordT>* >(valueRange);
	Enum* enumeration =  new Enum;
	for (typename std::set<WordT>::iterator it = vRE->values.begin(); it != vRE->values.end(); ++it)
	{
		enumeration->insert(*it);
	}

	return RangePtr(enumeration);
}


} // namespace Io
} // namespace CCfg

#endif /* CCfg_Io_BinaryFormat_h_INCLUDED */
