// $Id$
#ifndef CCFG_IO_DB_RUNCONDITION_H_INCLUDED
#define CCFG_IO_DB_RUNCONDITION_H_INCLUDED

#include "ConfigEntity.h"
#include "EntitySets.h"
#include "CCfg/Attribute.h"
#include <boost/intrusive_ptr.hpp>
#include <memory>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// forward declarations
////////////////////////////////////////////////////////////////////////////////
struct RunCondition;
typedef boost::intrusive_ptr<RunCondition> RunConditionPtr;
struct Run;
typedef boost::intrusive_ptr<Run> RunPtr;
struct ConfigVersion;
typedef boost::intrusive_ptr<ConfigVersion> ConfigVersionPtr;
////////////////////////////////////////////////////////////////////////////////
// class declarations
////////////////////////////////////////////////////////////////////////////////
/**
 * This entity stands for a reference to a configuration for a run.
 *
 * @author psizun
 */
struct RunCondition: public ConfigEntity
{
	static const char TABLE_NAME[];
	enum ColumnIndexes {ID = 0, FK_RUN, FK_CONFIGVERSION};

	UID fkRun;
	UID fkConfigVersion;

	RunCondition();
	virtual ~RunCondition();

	UID getId() const { return id; }

	bool operator==(const RunCondition& entity);
	bool operator!=(const RunCondition& entity) { return not operator==(entity); }

	RunPtr getRun();
	ConfigVersionPtr getConfigVersion();

private:
	RunCondition(const RunCondition& rvalue);
	RunCondition& operator=(const RunCondition& rvalue);
};

typedef boost::multi_index_container<
	RunConditionPtr,
	indexed_by<
		ordered_unique<tag<Id>, BOOST_MULTI_INDEX_CONST_MEM_FUN(RunCondition, UID, getId)>
	>
> RunConditionPtrSet;
typedef std::auto_ptr<RunConditionPtrSet>	RunConditionPtrSetPtr;
typedef RunConditionPtrSet::index<Id>::type RunConditionPtrIdIdx;

}
}
}

#endif /* CCFG_IO_DB_RUNCONDITION_H_INCLUDED */
