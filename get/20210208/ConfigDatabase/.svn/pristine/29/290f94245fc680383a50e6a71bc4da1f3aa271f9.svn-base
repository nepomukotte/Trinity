#include "InsertionProcess.h"
#include "AbstractUpdateController.h"
#include "ConfigClass.h"
#include "CCfg/ValidationVisitor.h"
#include "CCfg/Attribute.h"
#include <stdexcept>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::InsertionProcess class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor. Initializes the insertion parameters.
 */
InsertionProcess::InsertionProcess(Config* cfgEntity,
	Attribute* configRoot,
	AbstractUpdateController& controller)
	: cfgEntity(cfgEntity), configRoot(configRoot), controller(controller),
	pkAllocator(cfgEntity->getDAO()),
	newClassesBatch(), newClassInfosBatch(), newAttributesBatch(),
	classes(new ConfigClassPtrSet), classInfos(new ConfigClassInfoPtrSet)
{
}

/**
 * Destructor.
 */
InsertionProcess::~InsertionProcess()
{
}

/**
 * Executes the insertion process.
 * @return Returns false if the user cancelled the structure update.
 */
bool InsertionProcess::insert(bool skipValidation)
{
	if(not skipValidation)
	{
		ValidationVisitor validator;
		configRoot->accept(validator);
	}
	// retrieves the current structure version.
	ConfigStructureVersionPtr currentVersion = cfgEntity->getLastStructureVersion();
	if(currentVersion.get() != 0)
	{
		// a structure already exists: we have to prepare a structure update
		if(not prepareStructureUpdate(currentVersion.get()))
		{
			return false;
		}
	}
	else
	{
		// the configuration is empty: we have to create its first structure version.
		prepareFirstStructure();
	}
	// structure version is prepared (lists classes and class infos are filled)
	// we must create the new version.
	newVersion = structVersion->createNewVersion();
	newVersionsBatch.append(*newVersion);

	fillAttributesLists();
	insertPreparedData();
	return true;
}

/**
 * Prepares a structure update by comparing the new structure with the old one,
 * and by filling classes and classes info accordingly.
 */
bool InsertionProcess::prepareStructureUpdate(ConfigStructureVersion* oldStructure)
{
	// Retrieves classes from the configuration. New classes will be mixed with
	// these ones to have the full list of classes for the current configuration.
	// We not only retrieve classes from the previous structure, but also from
	// the whole configuration because we want to avoid creating a new distinct
	// class for a restored class.
	classes = cfgEntity->getClasses();
	// retrieves classes informations. This is necessary to reconstruct the old structure.
	ConfigClassInfoPtrSetPtr oldClassInfos = oldStructure->getClassInfos();

	std::auto_ptr<StructureNode> newStruct = StructureNode::generateStructure(*configRoot);
	std::auto_ptr<StructureNode> oldStruct = oldStructure->buildStructure(*classes, *oldClassInfos);

	StructureComparator comparator(classes, controller);
	bool isEqual = StructureNode::compare(*newStruct, *oldStruct, comparator);
	if(not isEqual)
	{
		if(controller.doConfirmStructureChange())
		{
			// structure version changed, we have to create a new record for it
			structVersion = cfgEntity->createNewStructure();
			newStructuresBatch.append(*structVersion);
			// make the list of new classes in comparison with the old version.
			fillNewClassesLists(comparator.getNewClasses());
			// recreate the class info list
			fillClassInfosLists(*newStruct);
		}
		else // the controller refused the structure update.
		{
			return false;
		}
	}
	else // structure have not changed
	{
		structVersion = oldStructure;
		// keep the class list and the class info list as it is.
		classInfos = oldClassInfos;
	}
	return true;
}

/**
 * Fills the class lists and the class information lists for a first structure insertion.
 */
void InsertionProcess::prepareFirstStructure()
{
	std::auto_ptr<StructureNode> newStruct = StructureNode::generateStructure(*configRoot);
	structVersion = cfgEntity->createNewStructure();
	newStructuresBatch.append(*structVersion);
	fillClassesLists(*newStruct);
	fillClassInfosLists(*newStruct);
}

/**
 * Fills the classes lists based on the result of the comparison between the current
 * and the new structure.
 */
void InsertionProcess::fillNewClassesLists(const NewClassList& newClasses)
{
	PrimaryKeyGenerator::SingleRange range;
	pkAllocator.reserve(range, TClass, newClasses.size());
	for(size_t i = 0; i < newClasses.size(); ++i)
	{
		const StructureNode& newClass = *newClasses[i];
		ConfigClassPtr cls(new ConfigClass);
		cls->id = range.nextVal();
		cls->fkCreationStructure = structVersion->id;
		cls->name = newClass.name;
		cls->path = newClass.path;

		if(newClass.parent == 0)
		{
			cls->fkParent = NULL_UID;
		}
		else
		{
			ConfigClassPathIdx& pathIdx = classes->get<Path>();
			ConfigClassPathIdx::const_iterator itParent =
				pathIdx.find(newClass.parent->path);
			if(itParent == pathIdx.end())
			{
				throw std::logic_error("The config class \"" +
						newClass.parent->path + "\" is missing.");
			}
			cls->fkParent = (*itParent)->id;
		}

		classes->insert(cls);
		newClassesBatch.append(*cls);
	}
}

/**
 * Creates the list of classes from a structure tree.
 */
void InsertionProcess::fillClassesLists(const StructureNode& structure)
{
	PrimaryKeyGenerator::SingleRange range;
	pkAllocator.reserve(range, TClass, StructureNode::computeNodeCount(structure));
	fillClassesLists(range, structure, NULL_UID);
}

/**
 * Recursively fills the list of classes from a structure tree.
 */
void InsertionProcess::fillClassesLists(PrimaryKeyGenerator::SingleRange& range,
	const StructureNode& currentClass,
	UID parentId)
{
	ConfigClassPtr cls(new ConfigClass);
	cls->id = range.nextVal();
	cls->fkParent = parentId;
	cls->fkCreationStructure = structVersion->id;
	cls->name = currentClass.name;
	cls->path = currentClass.path;

	classes->insert(cls);
	newClassesBatch.append(*cls);

	for(StructureNode::List::const_iterator it = currentClass.children.begin(),
		end = currentClass.children.end(); it != end; ++it)
	{
		fillClassesLists(range, *it, cls->id);
	}
}

/**
 * Creates the list of class informations from a structure tree.
 */
void InsertionProcess::fillClassInfosLists(const StructureNode& newStructure)
{
	PrimaryKeyGenerator::SingleRange range;
	pkAllocator.reserve(range, TClassInfo, StructureNode::computeNodeCount(newStructure));
	fillClassInfosLists(range, newStructure);
}

/**
 * Recursively fills the list of class informations from a structure tree.
 */
void InsertionProcess::fillClassInfosLists(PrimaryKeyGenerator::SingleRange& range,
	const StructureNode& currentClass)
{
	ConfigClassInfoPtr info(new ConfigClassInfo);
	info->id = range.nextVal();
	info->fkStructure = structVersion->id;
	info->type = currentClass.type;
	info->unit = currentClass.unit;
	info->range = currentClass.range;
	// retrieves the class associated with the current path.
	ConfigClassPathIdx& pathIdx = classes->get<Path>();
	ConfigClassPathIdx::const_iterator clsIt = pathIdx.find(currentClass.path);
	if(clsIt == pathIdx.end())
	{
		throw std::logic_error("The config class \"" +
				currentClass.path + "\" is missing.");
	}
	info->fkClass = (*clsIt)->id;
	classInfos->insert(info);
	newClassInfosBatch.append(*info);

	for(StructureNode::List::const_iterator it = currentClass.children.begin(),
		end = currentClass.children.end(); it != end; ++it)
	{
		fillClassInfosLists(range, *it);
	}
}

/**
 * Creates the list of attributes.
 */
void InsertionProcess::fillAttributesLists()
{
	PrimaryKeyGenerator::SingleRange range;
	pkAllocator.reserve(range, TAttribute, Attribute::computeAttributeCount(*configRoot));
	fillAttributesLists(range, *configRoot, NULL_UID);
}

/**
 * Recursively fills the list of attributes.
 */
void InsertionProcess::fillAttributesLists(PrimaryKeyGenerator::SingleRange& range,
	const Attribute& currentAttr, UID parentId)
{
	ConfigAttributePtr attr(new ConfigAttribute);
	attr->id = range.nextVal();
	attr->fkParent = parentId;
	attr->fkVersion = newVersion->id;
	attr->index = currentAttr.getIndex();
	if(not currentAttr.isObject())
	{
		attr->value = currentAttr.getValueAsString();
	}

	std::string path = currentAttr.getShortClassPath().toString();
	ConfigClassPathIdx& pathIdx = classes->get<Path>();
	ConfigClassPathIdx::const_iterator itClass = pathIdx.find(path);
	if(itClass == pathIdx.end())
	{
		throw std::logic_error("The config class \"" + path + "\" is missing.");
	}
	ConfigClassInfoClassIdx& classIdx = classInfos->get<Cfgclass>();
	ConfigClassInfoClassIdx::const_iterator itInfo = classIdx.find((**itClass).id);
	if(itInfo == classIdx.end())
	{
		throw std::logic_error("Missing class info associated with the class \"" + path + "\"");
	}

	attr->fkClassInfo = (**itInfo).id;
	newAttributesBatch.append(*attr);

	if(currentAttr.isObject())
	{
		for(Attribute::const_iterator itAttr = currentAttr.begin(),
			endAttr = currentAttr.end(); itAttr != endAttr; ++itAttr)
		{
			fillAttributesLists(range, *itAttr, attr->id);
		}
	}
}

/**
 * Inserts the content of all entity batches of this object, filled by the previous
 * method calls.
 */
void InsertionProcess::insertPreparedData()
{
	cfgEntity->getDAO()->insertFullConfig(newStructuresBatch, newVersionsBatch,
			newClassesBatch, newClassInfosBatch, newAttributesBatch);
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Config::StructureComparator class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
InsertionProcess::StructureComparator::StructureComparator(ConfigClassPtrSetPtr& classes, AbstractUpdateController& insController)
	: insController(insController), classes(classes)
{
}

/**
 * Destructor.
 */
InsertionProcess::StructureComparator::~StructureComparator()
{
}

/*
 * (see inherited doc)
 */
void InsertionProcess::StructureComparator::onLeftClassNotFound(const StructureNode& rightClass)
{
	insController.onClassRemoved(rightClass);
}

/*
 * (see inherited doc)
 */
void InsertionProcess::StructureComparator::onRightClassNotFound(const StructureNode& leftClass)
{
	recursiveAppendNewClass(leftClass);
	insController.onClassAdded(leftClass);
}

/*
 * (see inherited doc)
 */
void InsertionProcess::StructureComparator::onDifferentUnit(const StructureNode& leftClass, const StructureNode& rightClass)
{
	insController.onUnitChanged(leftClass, rightClass);
}

/*
 * (see inherited doc)
 */
void InsertionProcess::StructureComparator::onDifferentRange(const StructureNode& leftClass, const StructureNode& rightClass)
{
	insController.onRangeChanged(leftClass, rightClass);
}

/*
 * (see inherited doc)
 */
void InsertionProcess::StructureComparator::onDifferentType(const StructureNode& leftClass, const StructureNode& rightClass)
{
	insController.onTypeChanged(leftClass, rightClass);
}

/**
 * Appends the specified class and all its children (recursively) to the list of new classes.
 */
void InsertionProcess::StructureComparator::recursiveAppendNewClass(const StructureNode& root)
{
	// look for the class among the list of existing classes of the configuration
	// because the class may have existed before, and now be restored.
	ConfigClassPathIdx& pathIdx = classes->get<Path>();
	ConfigClassPathIdx::const_iterator itClass = pathIdx.find(root.path);
	if(itClass == pathIdx.end())
	{
		newClasses.push_back(&root);
	}

	for(StructureNode::List::const_iterator curChild = root.children.begin(),
		endChild = root.children.end(); curChild != endChild; ++curChild)
	{
		recursiveAppendNewClass(*curChild);
	}
}

}
}
}
