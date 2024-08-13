/*
 * RDocumentI.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#include "RDocumentI.h"
#include "RCConfigI.h"
#include ICE_HEADER(ObjectAdapter.h)

using namespace std;

namespace CCfg
{
/**
 * Load a config from a document
 */
RCConfigPrx	RDocumentI::getConfig(const Ice::Current& current) {
	RCConfigPrx proxy = RCConfigPrx::uncheckedCast(
		current.adapter->addWithUUID(
			new RCConfigI(doc->getConfig())
		)
	);
	return proxy;
}

/**
 * Save config to specified file
 */
void RDocumentI::saveToFile(const string& fileName, const Ice::Current&) {
	doc->saveToFile(fileName);
}

/*!
 * Save config to specified url
 */
void RDocumentI::saveToUrl(const string& url, const Ice::Current&) {
	doc->saveToUrl(url);
}

/*!
 * Save current config
 */
void RDocumentI::save(const Ice::Current&) {
	doc->save();
}


}
