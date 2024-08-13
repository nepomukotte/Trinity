/*
 * RDocumentI.h
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#ifndef CCfg_RDocumentI_h_INCLUDED
#define CCfg_RDocumentI_h_INCLUDED

#include "ConfigServer.h"
#include "DocumentFacade.h"

namespace CCfg
{

class RDocumentI : public RDocument
{
	DocumentFacadePtr doc;

public:
	RDocumentI(DocumentFacadePtr doc) : doc(doc) {}
	RCConfigPrx	getConfig(const Ice::Current&);
	void 		save(const Ice::Current&);
	void		saveToFile(const std::string& fileName, const Ice::Current&);
	void		saveToUrl(const std::string& url, const Ice::Current&);
};

}

#endif // CCfg_RDocumentI_h_INCLUDED
