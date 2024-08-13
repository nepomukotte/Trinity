/*
 * CConfigFacade.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#include "CConfigFacade.h"
#include "CatchExceptions.h"
#include <CCfg/Object.h>
#include <CCfg/Range.h>

using namespace std;

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfigFacade class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Create a new Config from a Path
 */
CConfigFacadePtr CConfigFacade::subConf(const string& path) {
	CConfigFacadePtr subconf;
	HANDLE_EXCEP(subconf.reset(new CConfigFacade(conf.appendPath(path))))
	return subconf;
}


/**
 * Set a new Path
 */
void	CConfigFacade::append(const string& path) {
	HANDLE_EXCEP(currpath.append(Path(path)))
}

/**
 * Test if "path" Node exist
 */
bool	CConfigFacade::testNode(const string& path) {
	bool ret = 0;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(path).exists())
	return ret;
}

/**
 * Create the Node "path"
 */
void   	CConfigFacade::createNode(const string& path) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(path).create(OBJECT_INST))
}

/**
 * Remove the Specified Node
 */
void   	CConfigFacade::removeNode(const string& path) {
	HANDLE_EXCEP(
			CCfg::Path	lastPath;

			if (conf.appendPath(currpath).appendPath(path).exists())
			{
				Attribute& myNode = *conf.getBaseAttribute(); 				//Fill attributes
				if (path != "")
					currpath = path;
				myNode.remove(currpath.toString());
				currpath.clear();
			}
		)
}

/**
 * Enumerate the Specified Node
 */
void 	CConfigFacade::enumNode(const string& path, AttributeList& nodeList) {
	HANDLE_EXCEP(
			CCfg::Path	lastPath;

			lastPath = currpath; 												//Save the current path to last path
			conf.appendPath(lastPath).appendPath(path);
			CConfig::Iterator it = conf.iterate();
			while(it.hasNext())
			{
				CConfig& childConf = it.next();
				PathNode node(childConf.find());
				if (childConf.getIndex() != "")
					nodeList.push_back(node.name+"["+childConf.getIndex()+"]"); //Save each Attribute
				else
					nodeList.push_back(node.name);
			}
			currpath.clear(); 					//Clear the Path
			currpath.append(lastPath); 			//Go Back to last Path
		)
}

/**
 * Get an Int Value
 */
int		CConfigFacade::getInt(const string& item) {
	int ret = 0;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getIntValue())
	return ret;
}

/**
 * Returns the integer value of an hexadecimal scalar attribute.
 * @param relativePath Relative path of the attribute.
 */
unsigned CConfigFacade::getHex(const string & item) {
	unsigned int ret = 0;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getHexValue())
	return ret;
}

/**
 * Returns the conversion to integer of the value in the concrete attribute.
 * @param relativePath Relative path of the attribute.
 */
int		CConfigFacade::getAsInt(const string & relativePath) {
	int ret = 0;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(relativePath).getValueAsInt())
	return ret;
}

/**
 * Returns the conversion to string of the value in the concrete attribute.
 * @param relativePath Relative path of the attribute.
 */
string	CConfigFacade::getAsStr(const string & relativePath) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(relativePath).getValueAsString())
	return ret;
}

/**
 * Get a String Value
 */
string	CConfigFacade::getStr(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getStringValue())
	return ret;
}

/**
 * Get a Bool Value
 */
bool 	CConfigFacade::getBool(const string& item) {
	bool ret = false;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getBoolValue())
	return ret;
}

/**
 * Get a Real Value
 */
float 	CConfigFacade::getReal(const string& item) {
	float ret = 0.0;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getRealValue())
	return ret;
}

/**
 * Get the Unit
 */
string	CConfigFacade::getUnit(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getUnit())
	return ret;
}

/**
 * Get the Name
 */
string	CConfigFacade::getName(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getName())
	return ret;
}

/**
 * Get the Index
 */
string	CConfigFacade::getIndex(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getIndex())
	return ret;
}

/**
 *
 */
string	CConfigFacade::getRange(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getRange()->toString())
	return ret;
}

/**
 * Get the Pre comment
 */
string	CConfigFacade::getPreComment(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getPreComment().getText())
	return ret;
}

/**
 * Get the Post comment
 */
string	CConfigFacade::getPostComment(const string& item) {
	string ret;
	HANDLE_EXCEP(ret = conf.appendPath(currpath).appendPath(item).getPostComment().getText())
	return ret;
}

/**
 * Get the current Path from the Root level
 */
string	CConfigFacade::getCurrentPath() {
	string ret;
	HANDLE_EXCEP(ret = currpath.toString())
	return ret;
}

/**
 * Sets an Int attribute from a string value.
 */
void 	CConfigFacade::setInt(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setIntValue(value))
}

/**
 * Sets an Int attribute from an integer value.
 */
void 	CConfigFacade::setInt(const string& item,  int value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setIntValue(value))
}

/**
 * Set a String value
 */
void 	CConfigFacade::setStr(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setStringValue(value))
}

/**
 * Sets a Bool attribute from a string value.
 */
void 	CConfigFacade::setBool(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setBoolValue(value))
}

/**
 * Sets a Bool attribute from a boolean value.
 */
void CConfigFacade::setBool(const string& item, bool value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setBoolValue(value))
}

/**
 * Sets a Real attribute from a string value.
 */
void CConfigFacade::setReal(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setRealValue(value))
}

/**
 * Sets a Real attribute from a floating-point value.
 */
void 	CConfigFacade::setReal(const string& item, float value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setRealValue(value))
}

/**
 * Sets an Hex attribute from its string value
 */
void 	CConfigFacade::setHex(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setHexValue(value))
}

/**
 * Sets an Hex value from its integer value
 */
void CConfigFacade::setHex(const string & item, unsigned value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setHexValue(value))
}

/**
 * Set the Unit
 */
void 	CConfigFacade::setUnit(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setUnit(value))
}

/**
 * Set the Unit
 */
void 	CConfigFacade::setName(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setName(value))
}

/**
 * Set of the Name
 */
void 	CConfigFacade::setIndex(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setIndex(value))
}

/**
 * Set of the Index.
 */
void 	CConfigFacade::setRange(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setRange(Range::createFromString(value)))
}

/**
 * Set the Pre comment
 */
void 	CConfigFacade::setPreComment(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setPreComment(value))
}

/**
 * Set the Post comment
 */
void 	CConfigFacade::setPostComment(const string& item, const string& value) {
	HANDLE_EXCEP(conf.appendPath(currpath).appendPath(item).setPostComment(value))
}

/**
 * Set path to "item" from root level
 */
void 	CConfigFacade::setCurrentPath(const string& item) {
	HANDLE_EXCEP(
			string st1 = ""; string st2 = "Undefined";
			currpath.clear();conf.appendPath(item);currpath = item;
			if (currpath == st1) { currpath = st2; }
		)
}


}
