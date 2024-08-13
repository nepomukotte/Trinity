// =============================================================================
// BinaryFormat.cpp
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

#include "BinaryFormat.h"
#include "../ValueRange.h"
#include "CCfg/Io/StreamStorage.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/TypeInfo.h"
#include "FileStorage.h"
#include "Url.h"
#include <boost/version.hpp>
#include <string>
#include <utility>
#include <fstream>
#include <sstream>
#include <stdexcept>

using std::pair;
using std::string;

namespace CCfg
{

namespace Io
{
/**
 * Constructor.
 */
BinaryFormat::BinaryFormat()
	: version(1), stringCount(0), attributeClassCount(0)
{
 cleanUp();
}

//*********************** FORMAT METHODS ***********************//

/**
 * Writes a configuration tree in a storage.
 * @param config the configuration tree.
 * @param storage the storage where data is written.
 * @param handler a callback to react to include directives.
 * @throws Storage::OpenError if there was an error opening the storage.
 */
void BinaryFormat::format(const Object* config, Storage& storage, IncludeHandler& /*handler*/)
{
	StreamStorage& stream = dynamic_cast<StreamStorage&>(storage);
	cleanUp();
	format(config);

	for(ByteBuffer::const_iterator it = ioBuffer.begin(),
			itEnd = ioBuffer.end(); it != itEnd; ++it)
	{
		stream.write(&*it, sizeof(Byte));
	}
	stream.close();
}
/**
 * The real format specialization called by the other format
 * @param attr Root attribute of the configuration tree
 */
void BinaryFormat::format(const Attribute* attr)
{
	recursiveFormat(attr);
	repack();
}

/**
 * Recursive format which moves in the CConfig Tree on each Attribute
 * @param attr Attribute processed at each loop
 */
void BinaryFormat::recursiveFormat(const Attribute* attr)
{
	processAttribute(attr); //process the current attribute
	if(attr->isObject())
	{
		for(Attribute::const_iterator it = attr->begin(); it != attr->end(); ++it)
		{
			recursiveFormat(it); //call the child attribute
		}
	}
}

/**
 * Transform an Attribute into an AttributeClass and process it.
 * Then it serialize the Attribute
 * @param attr Attribute to process
 */
void BinaryFormat::processAttribute(const Attribute* attr)
{
	AttributeClass attrClass(attr); // Creates corresponding Attribute class
	AttributeClassId classId = processAttributeClass(&attrClass);
	serializeAttribute(attr, classId);
}
/**
 * Put an Attribute Class into the AttributeClassToIdMap and serialize it if necessary
 * @param attrClass Attribute Class to process
 * @return Id corresponding to its place in the AttributeClassToIdMap
 */
AttributeClassId BinaryFormat::processAttributeClass(const AttributeClass* attrClass)
{
	pair<AttributeClassToIdMap::iterator, bool> attrClassRef;
	attrClassRef = attributeClassToIdMap.insert (
		make_pair( attrClass->getPath(), static_cast<AttributeClassId>(attributeClassToIdMap.size()) )
	);

	if (attrClassRef.second) // If AttributeClass not inserted in table yet
	{
		serializeAttributeClass(attrClass);
		attributeClassCount++;
	}

	return attrClassRef.first->second;
}

/**
 * Serialize the AttributeClass in the attributeClassBuffer.
 * Attribute only does not have all the info needed
 * @param attrClass Attribute Class to serialize
 */
void BinaryFormat::serializeAttributeClass(const AttributeClass* attrClass)
{
	//Insert the Class name
	StringId strId = processString(attrClass->getName());
	serializePrimitive<StringId, StringId>(attributeClassBuffer, strId);

	//Type insertion
	serializePrimitive<AttributeTypeWord, AttributeTypeId>(attributeClassBuffer, attrClass->getTypeId());

	//Unit insertion
	strId = processString(attrClass->getUnit());
	serializePrimitive<StringId, StringId>(attributeClassBuffer, strId);

	//Range insertion
	serializeRange(attributeClassBuffer, attrClass->getValueRange());
}

/**
 * Range serialization, 3 cases:
 * - No range
 * - Interval
 * - Enum
 * @param buff where the range is serialized
 * @param vRange Range to serialize
 *
 */
void BinaryFormat::serializeRange(ByteBuffer& buff, const ValueRange* vRange)
{
	RangeCodeWord rcw = (vRange == 0) ? NO_RANGE : vRange->code;
	serializePrimitive<RangeCodeWord, RangeCodeWord>(buff, rcw); //NO_RANGE, INTERVAL, ENUM

	switch(rcw)
	{
	//Interval template specialization
	case INTERVAL:
		switch (vRange->valueType)
		{
		case Id<Int>::val:
			serializeBoundaries<Int,IntWord>(vRange);
			break;

		case Id<Hex>::val:
			serializeBoundaries<Hex,HexWord>(vRange);
			break;

		case Id<Real>::val:
			serializeBoundaries<Real,RealWord>(vRange);
			break;

		case Id<Bool>::val:
			serializeBoundaries<Bool,BoolWord>(vRange);
			break;

		case Id<String>::val:
			serializeBoundaries<String,StringId>(vRange);
			break;

		default:
			break;
		}
		break;

		//Enum template specialization
	case ENUM:
		switch (vRange->valueType)
		{
		case Id<Int>::val:
			serializeEnum<Int,IntWord>(vRange);
			break;

		case Id<Hex>::val:
			serializeEnum<Hex,HexWord>(vRange);
			break;

		case Id<Real>::val:
			serializeEnum<Real,RealWord>(vRange);
			break;

		case Id<Bool>::val:
			serializeEnum<Bool,BoolWord>(vRange);
			break;

		case Id<String>::val:
			serializeEnum<String,StringId>(vRange);
			break;

		default:
			break;
		}
		break;
	default:
		break;
	}
}
/**
 * Serialize an Attribute:
 * - ClassId
 * - Attribute Index
 * - If it is an Object, it writes its ChildrenCount else its value
 * @param attr Attribute to serialize
 * @attrClassId Id associated to the Attribute
 */
void BinaryFormat::serializeAttribute(const Attribute* attr, AttributeClassId attrClassId)
{
	//Class Id serialization
	serializePrimitive<AttributeClassId,AttributeClassId>(attributeBuffer,attrClassId);
	//Attribute Index serialization
	serializePrimitive<StringId,String>(attributeBuffer,attr->getIndex());

	if(attr->getType() == OBJECT) //serialize child count
	{
		ChildCount childCount = 0;
		for(Attribute::const_iterator it = attr->begin(); it != attr->end(); it ++)
			childCount ++;

		serializePrimitive<ChildCount,ChildCount>(attributeBuffer,childCount);
	}
	else if(attr->getType() == INT_SCALAR)
	{
		serializePrimitive<IntWord,Int>(attributeBuffer,attr->getIntValue());
	}
	else if(attr->getType() == HEX_SCALAR)
	{
		serializePrimitive<HexWord,Hex>(attributeBuffer,attr->getHexValue());
	}
	else if(attr->getType() == BOOL_SCALAR)
	{
		serializePrimitive<BoolWord,Bool>(attributeBuffer,attr->getBoolValue());
	}
	else if(attr->getType() == REAL_SCALAR)
	{
		serializePrimitive<RealWord,Real>(attributeBuffer,attr->getRealValue());
	}
	else if(attr->getType() == STRING_SCALAR)
	{
		StringId id = processString(attr->getStringValue());
		serializePrimitive<StringId,StringId>(attributeBuffer,id);
	}
}

/**
 * Insert s string in the String map if necessary
 * Serialize the stringId in the StringTableBuffer
 * @param str String to process
 * @return Id relative to the string
 */
StringId BinaryFormat::processString(const std::string& str)
{
	// Add attribute name to string map (if not already done)
	pair<StringToIdMap::iterator, bool> stringRef;
	stringRef = stringToIdMap.insert (
		make_pair( str, static_cast<StringId>(stringToIdMap.size()) )
	);

	//If not in the table yet
	if (stringRef.second)
	{
		serialize2cstring(stringTableBuffer, str);
		stringCount++;
	}
	return stringRef.first->second;
}

/**
 * Serializes a string into a byte buffer in a form of
 * a null terminated sequence of characters.
 * @param buff where the string is transformed
 * @param str string to transform
 */
void BinaryFormat::serialize2cstring(ByteBuffer& buff, const std::string& str)
{
	for(size_t i=0; i<str.size(); i++)
	{
		buff.push_back(static_cast<Byte>(str[i]));
	}
	buff.push_back('\0');
}

/**
 * Function used to reassemble all the buffers in ioBuffer
 */
void BinaryFormat::repack()
{
	ioBuffer.clear();
	//Add the version number
	serializePrimitive<VersionNumber,VersionNumber>(ioBuffer,version);

	//Add the number of string
	serializePrimitive<TableSize,TableSize>(ioBuffer,stringCount);

	//Add the string table buffer
	copy(stringTableBuffer.begin(), stringTableBuffer.end(), back_inserter(ioBuffer));

	//Add the number of Attribute Class
	serializePrimitive<TableSize,TableSize>(ioBuffer,attributeClassCount);

	//Add the attribute class buffer
	copy(attributeClassBuffer.begin(), attributeClassBuffer.end(), back_inserter(ioBuffer));

	//Add the attribute buffer
	copy(attributeBuffer.begin(), attributeBuffer.end(), back_inserter(ioBuffer));
}

/**
 * Getter of config name
 * @return  The config name
 */
std::string BinaryFormat::getConfigName(Storage& storage)
{
	FileStorage* file = dynamic_cast<FileStorage*>(&storage);
	if(file != 0)
	{
		const std::string& path = file->getFileName();
		size_t pos = path.find_last_of(Url::PATH_SEPARATOR);
		if(pos != std::string::npos)
		{
			return path.substr(pos + 1);	// return the substring of the file name.
		}
		else // probably only a file name in the current directory.
		{
			return path;
		}
	}
	return "unnamed";
}

/**
 * Setter
 */
void BinaryFormat::setVersion(VersionNumber n)
{
	version = n;
}
/**
 * Cleaner. Erase all the previous value in all the buffers and maps
 */
void BinaryFormat::cleanUp()
{
	//Byte Buffer cleanup
	stringTableBuffer.clear();
	attributeClassBuffer.clear();
	attributeBuffer.clear();
	ioBuffer.clear();

	//variables cleanup
	stringCount = 0;
	attributeClassCount = 0;

	//Maps cleanup
	attributeClassToIdMap.clear();
	stringToIdMap.clear();
}

//*********************** FILE READER METHODS ***********************//

/**
 * Read a file from a string and fill the ioBuffer
 * @param url the file url
 */
void BinaryFormat::readBinaryFromFile(std::string url)
{
	cleanUp(); //cleans before starting
	int length = 0;
	char* buffer = 0;

	std::ifstream is(url.c_str(), std::ios::in | std::ios::binary );

	// get length of file:
	is.seekg (0, std::ios::end);
	length = is.tellg();
	is.seekg (0, std::ios::beg);

	// allocate memory:
	buffer = new char[length];
	// read data as a block:
	is.read (buffer,length);
	is.close();

	// Filling the ioBuffer
	for(int i=0 ; i<length ; ++i)
	{
		ioBuffer.push_back(buffer[i]);
	}

	delete[] buffer;

}

//*********************** PARSE METHODS ***********************//

/**
 * Parses the storage data to produce a configuration tree.
 * @param storage the storage where is the data.
 * @param handler a callback to react to include directives.
 * @return Returns the created configuration tree.
 * @throws Storage::OpenError if there was an error opening the storage.
 */
std::auto_ptr<CCfg::Object> BinaryFormat::parse(Storage& storage, IncludeHandler& /*handler*/)
{
	// Open storage
	StreamStorage *stream = dynamic_cast<StreamStorage*>(&storage);
	assert(stream != 0);

	if(not stream->isOpened())
	{
		stream->open(false);
	}

	// Read all data
	std::string data;
	stream->readAll(data);

	// Fill the buffer
	cleanUp();
	std::string::iterator byte;
	ioBuffer.resize(data.size());
	std::copy(data.begin(), data.end(), ioBuffer.begin());

	Object *attr = 0;
	try {
		//Attr is the main node of the config
		 attr = parse();
	}
	catch (std::out_of_range& e)
	{
		throw SemanticError();
	}

	return std::auto_ptr<CCfg::Object>(attr) ;
}

/**
 * Second implementation of parse, simpler version
 * @return Returns a pointer of the configuration tree
 */
Object* BinaryFormat::parse()
{

	size_t currentPos =0; //position in the ioBuffer
	version = parsePrimitive<VersionNumber,VersionNumber>(ioBuffer,currentPos);

	//Gets the number of strings
	TableSize stringTableNum = parsePrimitive<TableSize,TableSize>(ioBuffer,currentPos);

	//parse all the string of the ioBuffer
	parseStringMap(stringTableNum,currentPos);

	//Gets the number of Attribute class
	TableSize attributeClassNum = parsePrimitive<TableSize,TableSize>(ioBuffer,currentPos);

	//parse all the Attribute Class of the ioBuffer
	parseAttributeClassMap(attributeClassNum,currentPos);

	return dynamic_cast<Object*>(parseAttribute(ioBuffer,currentPos));
}

/**
 * Add the strings of the ioBuffer in the IdToStringMap
 * @param stringNum The number of string in the buffer
 * @param currentPos Current position in the buffers
 */
void BinaryFormat::parseStringMap(TableSize stringNum,size_t& currentPos )
{
	size_t count = 0;
	pair<IdToStringMap::iterator, bool> stringRef;

	while(count != stringNum)
		{
		//String created directly from the ByteBuffer
		std::string str(reinterpret_cast<const char*>(&ioBuffer[currentPos]));
		currentPos += str.size()+1; //updating the position in the buffer
		//insertion in the map
		stringRef = idToStringMap.insert (
				make_pair(static_cast<StringId>(idToStringMap.size()),str));
		count++;
		}
}

/**
 * Adds the new AttributeClass in the IdToAttributeClassMap
 * @param attrClassNum The number of attribute class in the buffer
 * @param currentPos Current position in the buffers
 */
void BinaryFormat::parseAttributeClassMap(TableSize attrClassNum,size_t& currentPos)
{
	for(size_t count = 0; count < attrClassNum; ++count)
	{
		std::auto_ptr<AttributeClass> attr = parseAttributeClass(currentPos); //gets an AttributeClass
		AttributeClassId nextId = static_cast<AttributeClassId>(idToAttributeClassMap.size());
		idToAttributeClassMap.insert(nextId, attr.release());
	}
}

/**
 * Rebuild an AttributeClass from the ioBuffer
 * @param buff where it extracts the Bytecode
 * @param currentPos Current position in the buffers
 * @return Pointer on the just parsed AttributeClass
 */
std::auto_ptr<AttributeClass> BinaryFormat::parseAttributeClass(size_t& currentPos)
{
	std::auto_ptr<AttributeClass> attr(new AttributeClass);
	AttributeTypeId attrType;
	RangeCodeWord rWord;
	ValueRangePtr rangePtr;

	//gets the name
	attr->setName(parsePrimitive<String,StringId>(ioBuffer,currentPos));

	 //the Attribute Type Word from the stream is converted in an AttributeTypeId
	attrType = parsePrimitive<AttributeTypeId,AttributeTypeWord>(ioBuffer,currentPos);
	attr->setTypeId(attrType);

	//Gets the unit
	attr->setUnit(parsePrimitive<String,StringId>(ioBuffer,currentPos));

	//Range parsing
	rWord=parsePrimitive<RangeCodeWord,RangeCodeWord>(ioBuffer,currentPos);
	if(rWord == NO_RANGE)
	{
		// Nothing
	}
	else if(rWord == INTERVAL)
	{
		switch(attrType)//Object, Int Real Hex Bool String
		{
		case INT_SCALAR:
			rangePtr.reset(parseRangeInterval<Int,IntWord>(ioBuffer,currentPos));
			break;

		case HEX_SCALAR:
			rangePtr.reset(parseRangeInterval<Hex,HexWord>(ioBuffer,currentPos));
			break;

		case REAL_SCALAR:
			rangePtr.reset(parseRangeInterval<Real,RealWord>(ioBuffer,currentPos));
			break;

		case BOOL_SCALAR:
			rangePtr.reset(parseRangeInterval<Bool,BoolWord>(ioBuffer,currentPos));
			break;

		case STRING_SCALAR:
			rangePtr.reset(parseRangeInterval<String,StringId>(ioBuffer,currentPos));
			break;

		default:
			break;
		}
	}
	else
	{
		switch(attrType)//Object, Int Real Hex Bool String
		{
		case INT_SCALAR:
			rangePtr.reset(parseRangeEnum<Int,IntWord>(ioBuffer,currentPos));
			break;

		case HEX_SCALAR:
			rangePtr.reset(parseRangeEnum<Hex,HexWord>(ioBuffer,currentPos));
			break;

		case REAL_SCALAR:
			rangePtr.reset(parseRangeEnum<Real,RealWord>(ioBuffer,currentPos));
			break;

		case BOOL_SCALAR:
			rangePtr.reset(parseRangeEnum<Bool,BoolWord>(ioBuffer,currentPos));
			break;

		case STRING_SCALAR:
			rangePtr.reset(parseRangeEnum<String,StringId>(ioBuffer,currentPos));
			break;

		default:
			break;
		}
	}

	attr->setValueRange(rangePtr);
	return attr;
}

/**
 * Rebuild an Attribute from a bytebuffer and the attribute class
 * @param buff where the data are parsed from
 * @param currentPos current position in the buffer
 * @return pointer on the freshly created attribute
 */
Attribute* BinaryFormat::parseAttribute(ByteBuffer& buff,size_t& currentPos)
{
	AttributeClass *attrClass;
	AttributeClassId id;
	IdToAttributeClassMap::iterator itAttrClass;
	String index;
	//gets the id of the attributeClass
	id = parsePrimitive<AttributeClassId,AttributeClassId>(buff,currentPos);
	//finds the attribute class corresponding to the id
	itAttrClass = idToAttributeClassMap.find(id);
	if (itAttrClass == idToAttributeClassMap.end()) throw SemanticError();
	//updating the pointer on the real AttributeClass
	attrClass = itAttrClass->second;

	//gets the index string
	index = parsePrimitive<String, StringId>(buff,currentPos);

	AttributePtr attr;
	//The attribute is an object
	if(attrClass->getTypeId() == OBJECT)
	{
		attr.reset(new Object);
		ChildCount childCount = parsePrimitive<ChildCount,ChildCount>(buff,currentPos);

		for(ChildCount i = 0; i < childCount; ++i )
		{
			attr->append(parseAttribute(buff,currentPos));
		}
	}
	//The attribute is a Scalar
	else
	{
		attr.reset(new Scalar);

		if (attrClass->getUnit() != "")
		{
			attr->setUnit(attrClass->getUnit());
		}

		if(attrClass->getValueRange() != 0)
		{
			attr->setRange(reaRange(attrClass->getValueRange()));
		}

		switch(attrClass->getTypeId())//Object, Int Real Hex Bool StringId
		{
		case INT_SCALAR:
			attr->setIntValue(parsePrimitive<Int,IntWord>(buff,currentPos));
			break;

		case HEX_SCALAR:
			attr->setHexValue(parsePrimitive<Hex,HexWord>(buff,currentPos));
			break;

		case REAL_SCALAR:
			attr->setRealValue(parsePrimitive<Real,RealWord>(buff,currentPos));
			break;

		case BOOL_SCALAR:
			attr->setBoolValue(parsePrimitive<Bool,BoolWord>(buff,currentPos));
			break;

		case STRING_SCALAR:
			attr->setStringValue(parsePrimitive<String,StringId>(buff,currentPos));
			break;

		default:
			break;

		}
	}
	attr->setName(attrClass->getName());
	attr->setIndex(index);

	return attr.release();
}

/**
 * Transform a ValueRange into a Range
 * @param valueRange the valuerange to transform
 * @return pointer on the Range created
 */
RangePtr BinaryFormat::reaRange(ValueRange* valueRange)
{
	RangePtr range;
	if(valueRange->code == INTERVAL)
	{
		switch(valueRange->valueType)
		{
		case INTEGER:
			{
				RangePtr temp = makeInterval<Int>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case REAL:
			{
				RangePtr temp = makeInterval<Real>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case HEXADECIMAL:
			{
				RangePtr temp = makeInterval<Hex>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case BOOLEAN:
			{
				RangePtr temp = makeInterval<Bool>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case STRING:
			{
				RangePtr temp = makeInterval<String>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case VALUETYPE_INVALID:
			break;

		default:
			break;
		}
	}
	else // If ENUM
	{
		switch(valueRange->valueType)
		{
		case INTEGER:
			{
				RangePtr temp = makeEnum<Int>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case REAL:
			{
				RangePtr temp = makeEnum<Real>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
				break;

		case HEXADECIMAL:
			{
				RangePtr temp = makeEnum<Hex>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case BOOLEAN:
			{
				RangePtr temp = makeEnum<Bool>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case STRING:
			{
				RangePtr temp = makeEnum<String>(valueRange);
				range = temp; // hack for VxWorks which does not support assignment from return!!
			}
			break;

		case VALUETYPE_INVALID:
			break;

		default:
			break;
		}
	}

	return range;

}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::BinaryFormat::SemanticError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
BinaryFormat::SemanticError::SemanticError()
{
}

/**
 * Destructor.
 */
BinaryFormat::SemanticError::~SemanticError() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string BinaryFormat::SemanticError::formatMessage() const
{
	std::ostringstream message;
	message << "The binary document is not a valid configuration.";
	return message.str();
}

}//namespace Io
}//namespace CCfg
