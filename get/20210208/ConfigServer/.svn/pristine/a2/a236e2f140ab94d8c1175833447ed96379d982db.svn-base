/*
 * ConfigServer.h
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#ifndef CCfg_ConfigServerI_h_INCLUDED
#define CCfg_ConfigServerI_h_INCLUDED

#include "ConfigServer.h"
#include <string>

namespace CCfg
{
class ConfigServerI	: public ConfigServer {
public:
	virtual RDocumentPrx	loadDoc(const std::string& path, const Ice::Current&);
	virtual RDocumentPrx	newDoc(const Ice::Current&);
};

}

#endif // CCfg_ConfigServerI_h_INCLUDED
