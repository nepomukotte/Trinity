#include "ConfigDao.h"
#include "ConfigStructureVersion.h"
#include "../../../Utils/Date/StdTmHelpers.h"
#include <map>
#include <utility>
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigStructureVersion class methods body
////////////////////////////////////////////////////////////////////////////////

const char ConfigStructureVersion::TABLE_NAME[] = "CCfg_ConfigStructureVersion";

/**
 * Constructor.
 */
ConfigStructureVersion::ConfigStructureVersion()
{
}

/**
 * Destructor.
 */
ConfigStructureVersion::~ConfigStructureVersion()
{
}

/**
 * Returns the date of creation of this structure version.
 */
const std::tm& ConfigStructureVersion::getCreationDate() const
{
	return creationDate;
}

/**
 * Return the date of creation of this structure version as a time_t (number of
 * seconds since the Epoch).
 */
time_t ConfigStructureVersion::getCreationDateAsTimeT()
{
	return mktime(&creationDate);
}

/**
 * Equality operator.
 */
bool ConfigStructureVersion::operator==(const ConfigStructureVersion& entity)
{
	return id == entity.id and
		creationDate == entity.creationDate and
		fkConfig == entity.fkConfig;
}

/**
 * @fn bool ConfigStructureVersion::operator!=(const ConfigStructureVersion& entity)
 * Unequality operator.
 */

/**
 * Returns the list of configuration classes in use for this structure version.
 */
ConfigClassPtrSetPtr ConfigStructureVersion::getClasses()
{
	assert(dao != 0);
	return dao->getClassesOfStructure(*this);
}

/**
 * Returns the list of configuration class informations for this structure version.
 */
ConfigClassInfoPtrSetPtr ConfigStructureVersion::getClassInfos()
{
	assert(dao != 0);
	return dao->getClassInfosOfStructure(*this);
}

/**
 * Retrieves the configuration associated with this version, from the database.
 */
ConfigPtr ConfigStructureVersion::getConfig()
{
	assert(dao != 0);
	return dao->getConfig(fkConfig);
}

/**
 * Creates a new version record for this structure.
 */
ConfigVersionPtr ConfigStructureVersion::createNewVersion()
{
	assert(dao != 0);
	ConfigVersionPtr version(new ConfigVersion);
	version->setDAO(dao);
	version->id = dao->allocIdRangeForVersions(1) + 1;
	version->fkStructure = this->id;
	time_t currentTime = time(0);
	version->creationDate = *gmtime(&currentTime);
	// name left to null
	return version;
}

/**
 * Builds the StructureNode tree for the specified list of classes and class infos
 * previously retrieved from the database.
 */
std::auto_ptr<StructureNode>
ConfigStructureVersion::buildStructure(ConfigClassPtrSet& classes,
	ConfigClassInfoPtrSet& infos)
{
	// FIXME: untested
	std::auto_ptr<StructureNode> rootSNode(new StructureNode); // the structure root node

	typedef std::map<UID, StructureNode*> Id2NodeMap;
	Id2NodeMap id2node;	// map to associate a class id with the corresponding StructureNode

	ConfigClassInfoIdIdx& infoIds = infos.get<Id>();
	ConfigClassIdIdx& clsIds = classes.get<Id>();

	// loop on every ConfigClassInfo
	for(ConfigClassInfoIdIdx::const_iterator it = infoIds.begin(), end = infoIds.end();
		it != end; ++it)
	{
		ConfigClassInfo& clsInfo = **it;

		// find the class to which the class info belongs
		ConfigClassIdIdx::const_iterator itCls = clsIds.find(clsInfo.fkClass);
		assert(itCls != clsIds.end());
		ConfigClass& cls = **itCls;

		if(cls.fkParent == NULL_UID) // initialize the root
		{
			initStructureNode(*rootSNode, cls, clsInfo);
			id2node.insert(std::make_pair(cls.id, rootSNode.get()));
		}
		else
		{
			// find the parent of the class
			ConfigClassIdIdx::const_iterator itParent = clsIds.find(cls.fkParent);
			assert(itParent != clsIds.end());
			ConfigClass& parentCls = **itParent;

			// find the structure node corresponding to the parent class
			Id2NodeMap::const_iterator itNode = id2node.find(parentCls.id);
			assert(itNode != id2node.end());
			StructureNode& parentNode = *itNode->second;

			// create a structure node from the current class & class info
			StructureNode node;
			initStructureNode(node, cls, clsInfo);
			// add it to the list of children of its parent
			parentNode.append(node);
			id2node.insert(std::make_pair(cls.id, &parentNode.children.back()));
		}
	}
	return rootSNode;
}

void ConfigStructureVersion::initStructureNode(StructureNode& node, ConfigClass& cls, ConfigClassInfo& clsInfo)
{
	node.name = cls.name;
	node.path = cls.path;
	node.type = clsInfo.type;
	node.unit = clsInfo.unit;
	node.range = clsInfo.range;
}

}
}
}
