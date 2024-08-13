#ifndef CCFG_IO_DATABASEFORMATFACTORY_H_INCLUDED
#define CCFG_IO_DATABASEFORMATFACTORY_H_INCLUDED

#include "CCfg/Io/AbstractFormatFactory.h"
#include "DatabaseFormat.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class creates DatabaseFormat objects from URLs.
 *
 * @author fchateau
 */
class DatabaseFormatFactory : public AbstractFormatFactory
{
	UpdateControllerPtr updateController;
public:
	DatabaseFormatFactory();
	~DatabaseFormatFactory();

	std::auto_ptr<Format> createFormat(const Url& url) const;

	const UpdateControllerPtr& getUpdateController() const { return updateController; }
	void setUpdateController(const UpdateControllerPtr& controller) { updateController = controller; }
};

}
}

#endif // CCFG_IO_DATABASEFORMATFACTORY_H_INCLUDED
