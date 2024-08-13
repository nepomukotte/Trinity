/*
 * CConfigFacade.h
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#ifndef CCfg_CConfigFacade_h_INCLUDED
#define CCfg_CConfigFacade_h_INCLUDED

#include <CCfg/CConfig.h>
#include <CCfg/Path.h>

namespace CCfg
{
//forward declaration to make the following typedef work.
class CConfigFacade;


//typedef std::auto_ptr<CConfigFacade> CConfigFacadePtr;
typedef boost::shared_ptr<CConfigFacade> CConfigFacadePtr;
typedef std::vector<std::string> AttributeList;

class CConfigFacade
{
	CCfg::CConfig	conf;
	CCfg::Path		currpath;

public:
	CConfigFacade() {}
	CConfigFacade(CCfg::CConfig conf1) : conf(conf1) {}
	~CConfigFacade() {}

	//CONFIG FUNCTIONS
	CConfigFacadePtr	subConf(const std::string& path);
	void		removeNode(const std::string& path);
	bool		testNode(const std::string& path);
	void		append(const std::string& path);
	void		createNode(const std::string& path);
	void		enumNode(const std::string& path, AttributeList& EnumList);

	//GET FUNCTIONS
	std::string 	getCurrentPath();
	int				getInt(const std::string& item);
	std::string 	getStr(const std::string& item);
	bool 			getBool(const std::string& item);
	float			getReal(const std::string& item);
	unsigned		getHex(const std::string & item);
	std::string 	getAsStr(const std::string& item);
	int				getAsInt(const std::string & item);
	std::string 	getUnit(const std::string& item);
	std::string 	getName(const std::string& item);
	std::string 	getIndex(const std::string& item);
	std::string 	getRange(const std::string& item);
	std::string 	getPreComment(const std::string& item);
	std::string 	getPostComment(const std::string& item);

	//SET FUNCTIONS
	void	setCurrentPath(const std::string& item);
	void	setInt(const std::string& item, const std::string& value);
	void	setInt(const std::string& item, int value);
	void	setStr(const std::string& item, const std::string& value);
	void	setBool(const std::string& item, const std::string& value);
	void	setBool(const std::string& item, bool value);
	void	setReal(const std::string& item, const std::string& value);
	void	setReal(const std::string& item, float value);
	void	setHex(const std::string & item, const std::string & value);
	void	setHex(const std::string & item, unsigned value);
	void	setUnit(const std::string& item, const std::string& value);
	void	setName(const std::string& item, const std::string& value);
	void	setIndex(const std::string& item, const std::string& value);
	void	setRange(const std::string& item, const std::string& value);
	void 	setPreComment(const std::string& item, const std::string& value);
	void 	setPostComment(const std::string& item, const std::string& value);
};

}

#endif // CCfg_CConfigFacade_h_INCLUDED
