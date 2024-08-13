/*
 * RCConfigI.h
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#ifndef CCfg_RCConfigI_h_INCLUDED
#define CCfg_RCConfigI_h_INCLUDED

#include "ConfigServer.h"
#include "CConfigFacade.h"

namespace CCfg
{

class RCConfigI : public RCConfig
{
	CConfigFacadePtr conf;

public:
	RCConfigI(CConfigFacadePtr facade) : conf(facade) {}

	//CONFIG FUNCTIONS
	RCConfigPrx	subConf(const std::string& path, const Ice::Current&);
	void		removeNode(const std::string& path, const Ice::Current&);
	bool		testNode(const std::string& path, const Ice::Current&);
	void		append(const std::string& path, const Ice::Current&);
	void		createNode(const std::string& path, const Ice::Current&);
	void		enumNode(const std::string& path, RAttributeList& EnumList, const Ice::Current&);

	//GET FUNCTIONS
	std::string 	getCurrentPath(const Ice::Current&);
	int				getInt(const std::string& item, const Ice::Current&);
	std::string 	getStr(const std::string& item, const Ice::Current&);
	bool 			getBool(const std::string& item, const Ice::Current&);
	float			getReal(const std::string& item, const Ice::Current&);
	std::string 	getUnit(const std::string& item, const Ice::Current&);
	std::string 	getName(const std::string& item, const Ice::Current&);
	std::string 	getIndex(const std::string& item, const Ice::Current&);
	std::string 	getRange(const std::string& item, const Ice::Current&);
	std::string 	getPreComment(const std::string& item, const Ice::Current&);
	std::string 	getPostComment(const std::string& item, const Ice::Current&);

	//SET FUNCTIONS
	void	setCurrentPath(const std::string& item, const Ice::Current&);
	void	setInt(const std::string& item, const std::string& value, const Ice::Current&);
	void	setStr(const std::string& item, const std::string& value, const Ice::Current&);
	void	setBool(const std::string& item, const std::string& value, const Ice::Current&);
	void	setReal(const std::string& item, const std::string& value, const Ice::Current&);
	void	setUnit(const std::string& item, const std::string& value, const Ice::Current&);
	void	setName(const std::string& item, const std::string& value, const Ice::Current&);
	void	setIndex(const std::string& item, const std::string& value, const Ice::Current&);
	void	setRange(const std::string& item, const std::string& value, const Ice::Current&);
	void 	setPreComment(const std::string& item, const std::string& value, const Ice::Current&);
	void 	setPostComment(const std::string& item, const std::string& value, const Ice::Current&);
};
}

#endif // CCfg_RCConfigI_h_INCLUDED
