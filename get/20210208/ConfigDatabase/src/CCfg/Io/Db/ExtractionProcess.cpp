#include "ExtractionProcess.h"
#include "ConfigStructureVersion.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/ScalarClassInfo.h"
#include "CCfg/Range.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ExtractionProcess class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
ExtractionProcess::ExtractionProcess(ConfigVersion& version)
	: version(version), result(),
	classes(), classInfos(), attributes(),
	clsIdIdx(0), infoIdIdx(0), attrIdIdx(0)
{
}

/**
 * Destructor.
 */
ExtractionProcess::~ExtractionProcess()
{
}

/**
 * Executes the extraction process.
 * @return Returns the extracted configuration.
 */
std::auto_ptr<Attribute>& ExtractionProcess::extract()
{
	loadData();
	for(ConfigAttributeIdIdx::iterator itAttr = attrIdIdx->begin(),
		endAttr = attrIdIdx->end(); itAttr != endAttr; ++itAttr)
	{
		ConfigAttribute& attr = **itAttr;
		ConfigClassInfo& info = getClassInfo(attr.fkClassInfo);
		ConfigClass& cls = getClass(info.fkClass);

		// instanciate the attribute according to its type.
		std::auto_ptr<Attribute> newAttrPtr;
		if(info.type == OBJECT)
		{
			newAttrPtr.reset(new Object(cls.name, attr.index));
		}
		else
		{
			newAttrPtr.reset(new Scalar(cls.name, attr.index));
			newAttrPtr->setValue(TypeInfo::lookup(info.type), attr.value);
		}
		Attribute* newAttr = newAttrPtr.get();
		// put the address of the instanciated attribute in the ConfigAttribute entity.
		// It will be used to quickly associate the parent ConfigAttribute with
		// the corresponding CCfg::Attribute.
		attr.userData = newAttr;

		if(attr.fkParent == NULL_UID)
		{
			result = newAttrPtr;
		}
		else
		{
			ConfigAttribute& dbParentAttr = getAttribute(attr.fkParent);
			// Retrieving the parent Attribute through the userData field.
			// we assume that the parent has already been processed because
			// ConfigAttributes are processed by increasing id, and foreign keys must
			// reference existing id (preceding the considered one).
			Attribute* parentAttr = reinterpret_cast<Attribute*>(dbParentAttr.userData);
			parentAttr->append(newAttrPtr.release());
		}
		if(not (info.unit.empty() and info.range.empty()))
		{
			ScalarClassInfo* classInfo = dynamic_cast<ScalarClassInfo*>(newAttr->getClassInfo().get());
			if(classInfo == 0)
			{
				newAttr->setUnit(info.unit);
				newAttr->setRange(Range::createFromString(info.range));
			}
		}
	}
	return result;
}

/**
 * Loads all data related to the selected configuration version from the database.
 */
void ExtractionProcess::loadData()
{
	ConfigStructureVersionPtr structure = version.getStructure();
	classes = structure->getClasses();
	classInfos = structure->getClassInfos();
	attributes = version.getAttributes();

	clsIdIdx = &classes->get<Id>();
	infoIdIdx = &classInfos->get<Id>();
	attrIdIdx = &attributes->get<Id>();
}

/**
 * Returns the loaded class with the specified id.
 * @throws ExtractionProcess::Error if it cannot be found.
 */
ConfigClass& ExtractionProcess::getClass(UID id)
{
	ConfigClassIdIdx::iterator itClass = clsIdIdx->find(id);
	if(itClass == clsIdIdx->end())
	{
		throw Error();
	}
	return **itClass;
}

/**
 * Returns the loaded class info with the specified id.
 * @throws ExtractionProcess::Error if it cannot be found.
 */
ConfigClassInfo& ExtractionProcess::getClassInfo(UID id)
{
	ConfigClassInfoIdIdx::iterator itClassInfo = infoIdIdx->find(id);
	if(itClassInfo == infoIdIdx->end())
	{
		throw Error();
	}
	return **itClassInfo;
}

/**
 * Returns the loaded attribute with the specified id.
 * @throws ExtractionProcess::Error if it cannot be found.
 */
ConfigAttribute& ExtractionProcess::getAttribute(UID id)
{
	ConfigAttributeIdIdx::iterator itAttr = attrIdIdx->find(id);
	if(itAttr == attrIdIdx->end())
	{
		throw Error();
	}
	return **itAttr;
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ExtractionProcess::Error class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ExtractionProcess::Error::Error() : Exception()
{
}

/**
 * Copy constructor.
 */
ExtractionProcess::Error::Error(const Error& rvalue)
	: Exception(rvalue)
{
}

/**
 * Assignment operator.
 */
ExtractionProcess::Error&
ExtractionProcess::Error::operator=(const Error& rvalue)
{
	Exception::operator=(rvalue);
	return *this;
}

/*
 * (see inherited doc)
 */
std::string
ExtractionProcess::Error::formatMessage() const
{
	return "Error loading a database configuration: incoherent data";
}

}
}
}
