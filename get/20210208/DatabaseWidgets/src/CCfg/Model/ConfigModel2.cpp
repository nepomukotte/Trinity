#include "ConfigModel2.h"
#include "CCfg/Model/UndoableConfigModelAccessor.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModel2 class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigModel2::ConfigModel2(QObject* parent) : ConfigModel(parent)
{
}

/**
 * Destructor.
 */
ConfigModel2::~ConfigModel2()
{
}

/**
 * Saves the configuration in its associated database.
 */
void ConfigModel2::commit()
{
	Io::Document* doc = accessor->getDocument();
	doc->exportSubConfig(*doc->getConfig(), Io::Url(dbUrl)); // use export to avoid changing the bound url.
}

/**
 * Saves the configuration in the database designated by the url and associate
 * it with the database.
 */
void ConfigModel2::commitInto(const Io::Url& url)
{
	Io::Document* doc = accessor->getDocument();
	doc->exportSubConfig(*doc->getConfig(), url);
	dbUrl = url.format();
}

/**
 * Changes the resource associated with the document.
 * @param url the url of the new resource we want to associate to the document.
 */
void ConfigModel2::setDocumentUrl(const std::string& url)
{
	accessor->getDocument()->setUrl(url);
}

}
}
