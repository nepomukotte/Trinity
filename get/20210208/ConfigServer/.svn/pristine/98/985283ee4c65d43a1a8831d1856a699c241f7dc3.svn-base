/*
 * DocumentFacade.h
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#ifndef CCfg_DocumentFacade_h_INCLUDED
#define CCfg_DocumentFacade_h_INCLUDED

#include "CConfigFacade.h"
#include "CCfg/Io/Document.h"

namespace CCfg
{
class DocumentFacade;

//typedef std::auto_ptr<DocumentFacade> DocumentFacadePtr;
typedef boost::shared_ptr<DocumentFacade> DocumentFacadePtr;

class DocumentFacade
{
	Io::Document	docu;
	std::string		docName;

public:
	DocumentFacade(const std::string& fileName = "Undefined") { docName = fileName; }
	~DocumentFacade() {}
	static DocumentFacadePtr newDoc();
	static DocumentFacadePtr loadDoc(const std::string& fileName);
	CConfigFacadePtr	getConfig();
	void 				save();
	void				saveToFile(const std::string& fileName);
	void				saveToUrl(const std::string& url);
};

}

#endif // CCfg_DocumentFacade_h_INCLUDED
