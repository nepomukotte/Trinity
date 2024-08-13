/*
 * DocumentFacade.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#include "DocumentFacade.h"
#include "CatchExceptions.h"
#include <CCfg/Io/Url.h>

using namespace std;

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::DocumentFacade class methods body
///////////////////////////////////////////////////////////////////////////////

DocumentFacadePtr DocumentFacade::newDoc()
{
	return DocumentFacadePtr(new DocumentFacade());
}

DocumentFacadePtr DocumentFacade::loadDoc(const std::string& fileName)
{
	return DocumentFacadePtr(new DocumentFacade(fileName));
}

/**
 * Load a config from a document
 */
CConfigFacadePtr	DocumentFacade::getConfig() {
	CConfigFacadePtr proxy;

	HANDLE_EXCEP(
		if (docName != "Undefined")
			docu.loadFromFile(docName);
		proxy.reset(new CConfigFacade(docu.getConfig()));
	)
	return proxy;
}

/**
 * Save config to specified file
 */
void		DocumentFacade::saveToFile(const string& fileName) {
	HANDLE_EXCEP(docu.saveToFile(fileName))
}

/**
 * Save config to specified url
 */
void		DocumentFacade::saveToUrl(const string& url) {
	HANDLE_EXCEP(docu.saveToUrl(Io::Url(url)))
}

/**
 * Save current config
 */
void		DocumentFacade::save() {
	HANDLE_EXCEP(docu.save())
}


}
