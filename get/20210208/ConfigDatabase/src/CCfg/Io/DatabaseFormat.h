#ifndef CCFG_IO_DATABASEFORMAT_H_INCLUDED
#define CCFG_IO_DATABASEFORMAT_H_INCLUDED

#include "Db/AbstractUpdateController.h"
#include "CCfg/Io/Format.h"
#include <boost/shared_ptr.hpp>

namespace CCfg
{
namespace Io
{
typedef boost::shared_ptr<Db::AbstractUpdateController> UpdateControllerPtr;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the implementation of the Format interface for the database mapping.
 *
 * @author fchateau
 */
class DatabaseFormat : public Format
{
	UpdateControllerPtr updateController;
public:
	class MappingError;

	enum Operation { Insertion, Extraction };

	// @{ construction
	DatabaseFormat();
	~DatabaseFormat();
	// @}
	// @{ accessors
	const UpdateControllerPtr& getUpdateController() const { return updateController; }
	void setUpdateController(const UpdateControllerPtr& controller) { updateController = controller; }
	// @}
	// @{ interface implementation
	std::auto_ptr<CCfg::Object> parse(Storage& storage, IncludeHandler& handler);
	void format(const Object *config, Storage& storage, IncludeHandler& handler);
	std::string getConfigName(Storage& storage);
	// @}
private:
	DatabaseFormat(const DatabaseFormat& rvalue);
	DatabaseFormat& operator=(const DatabaseFormat& rvalue);
};



/**
 * This exception reports an error during a database mapping operation.
 */
class DatabaseFormat::MappingError : public Exception
{
	Operation op;
	std::string msg;

public:
	MappingError(Operation op, const std::string& message);
	MappingError(const MappingError& rvalue);
	~MappingError() throw() {}
	MappingError& operator=(const MappingError& rvalue);

private:
	std::string formatMessage() const;
};

}
}

#endif // CCFG_IO_DATABASEFORMAT_H_INCLUDED
