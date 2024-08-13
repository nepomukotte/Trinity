#include "ConfigServerI.h"
#include "RDocumentI.h"
#include "DocumentFacade.h"
#include ICE_HEADER(ObjectAdapter.h)

using namespace std;

namespace CCfg {
///////////////////////////////////////////////////////////////////////////////
// CCfg::ConfigServerI class methods body
///////////////////////////////////////////////////////////////////////////////
/**
 * Load a new document
 */
RDocumentPrx	ConfigServerI::loadDoc(const string& fileName, const Ice::Current& current) {
	return RDocumentPrx::uncheckedCast(
		current.adapter->addWithUUID(
			new RDocumentI(DocumentFacade::loadDoc(fileName))
		)
	);
}
/**
 * Create a new document
 */
RDocumentPrx	ConfigServerI::newDoc(const Ice::Current& current) {
	return RDocumentPrx::uncheckedCast(
		current.adapter->addWithUUID(
			new RDocumentI(DocumentFacade::newDoc())
		)
	);
}
}
