#ifndef CCFG_IO_DB_INSERTIONPROCESS_H_INCLUDED
#define CCFG_IO_DB_INSERTIONPROCESS_H_INCLUDED

#include "ConfigStructureVersionBatch.h"
#include "ConfigVersionBatch.h"
#include "ConfigClassBatch.h"
#include "ConfigClassInfoBatch.h"
#include "ConfigAttributeBatch.h"
#include "PrimaryKeyGenerator.h"
#include "ConfigDao.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
// forward declarations
struct AbstractUpdateController;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class encapsulates the process of translating a compound config into a
 * database form, and performing the insertion.
 *
 * @author fchateau
 */
class InsertionProcess
{
	typedef std::vector<const StructureNode*> NewClassList;

	Config* cfgEntity;
	Attribute* configRoot;
	AbstractUpdateController& controller;

	PrimaryKeyGenerator pkAllocator;

	// @{ batches of entity to insert (for entities where the expected count is above 1)
	ConfigStructureVersionBatch newStructuresBatch;
	ConfigVersionBatch newVersionsBatch;
	ConfigClassBatch newClassesBatch;
	ConfigClassInfoBatch newClassInfosBatch;
	ConfigAttributeBatch newAttributesBatch;
	// @}

	ConfigStructureVersionPtr structVersion;
	ConfigVersionPtr newVersion;

	// @{ entity sets for entities where indexed search will have to be performed as
	// part of the algorithm
	ConfigClassPtrSetPtr classes;
	ConfigClassInfoPtrSetPtr classInfos;
	// @}

public:
	InsertionProcess(Config* cfgEntity, Attribute* configRoot, AbstractUpdateController& controller);
	virtual ~InsertionProcess();

	bool insert(bool skipValidation = false);
	ConfigVersionPtr getNewVersion() { return newVersion; }

private:
	bool prepareStructureUpdate(ConfigStructureVersion* oldStructure);
	void prepareFirstStructure();

	void fillNewClassesLists(const NewClassList& newClasses);
	void fillClassesLists(const StructureNode& structure);
	void fillClassesLists(PrimaryKeyGenerator::SingleRange& range, const StructureNode& currentClass, UID parentId);
	void fillClassInfosLists(const StructureNode& newStructure);
	void fillClassInfosLists(PrimaryKeyGenerator::SingleRange& range, const StructureNode& currentClass);
	void fillAttributesLists();
	void fillAttributesLists(PrimaryKeyGenerator::SingleRange& range, const Attribute& currentAttr, UID parentId);

	void insertPreparedData();

	InsertionProcess(const InsertionProcess& rvalue);
	InsertionProcess& operator=(const InsertionProcess& rvalue);

	class StructureComparator;
};

/**
 * This class implements the ComparisonMonitor to compute structural differences
 * between the memory configuration and the latest database configuration.
 */
class InsertionProcess::StructureComparator : public StructureNode::ComparisonMonitor
{
public:
	StructureComparator(ConfigClassPtrSetPtr& classes, AbstractUpdateController& insController);
	~StructureComparator();
	void onLeftClassNotFound(const StructureNode& rightClass);
	void onRightClassNotFound(const StructureNode& leftClass);
	void onDifferentUnit(const StructureNode& leftClass, const StructureNode& rightClass);
	void onDifferentRange(const StructureNode& leftClass, const StructureNode& rightClass);
	void onDifferentType(const StructureNode& leftClass, const StructureNode& rightClass);

	const NewClassList& getNewClasses() const { return newClasses; }

private:
	void recursiveAppendNewClass(const StructureNode& root);

	AbstractUpdateController& insController;
	NewClassList newClasses;
	ConfigClassPtrSetPtr& classes;
};

}
}
}

#endif // CCFG_IO_DB_INSERTIONPROCESS_H_INCLUDED
