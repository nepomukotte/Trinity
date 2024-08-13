/*
 * RCConfigI.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: rvloic
 */

#include "RCConfigI.h"
#include ICE_HEADER(ObjectAdapter.h)

using namespace std;

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::RCConfigI class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Create a new Config from a Path
 */
RCConfigPrx	RCConfigI::subConf(const string& path, const Ice::Current& current) {
	return RCConfigPrx::uncheckedCast(
		current.adapter->addWithUUID(
			new RCConfigI(conf->subConf(path))
		)
	);
}

/**
 * Set a new Path
 */
void	RCConfigI::append(const string& path, const Ice::Current&) {
	conf->append(path);
}

/**
 * Test if "path" Node exist
 */
bool	RCConfigI::testNode(const string& path, const Ice::Current&) {
	return conf->testNode(path);
}

/**
 * Create the Node "path"
 */
void   	RCConfigI::createNode(const string& path, const Ice::Current&) {
	conf->createNode(path);
}

/**
 * Remove the Specified Node
 */
void   	RCConfigI::removeNode(const string& path, const Ice::Current&) {
	conf->removeNode(path);
}

/**
 * Enumerate the Specified Node
 */
void 	RCConfigI::enumNode(const string& path, RAttributeList& nodeList, const Ice::Current&) {
	conf->enumNode(path, nodeList);
}

/**
 * Get an Int Value
 */
int		RCConfigI::getInt(const string& item, const Ice::Current&) {
	return conf->getInt(item);
}

/**
 * Get a String Value
 */
string	RCConfigI::getStr(const string& item, const Ice::Current&) {
	return conf->getStr(item);
}

/**
 * Get a Bool Value
 */
bool 	RCConfigI::getBool(const string& item, const Ice::Current&) {
	return conf->getBool(item);
}

/**
 * Get a Real Value
 */
float 	RCConfigI::getReal(const string& item, const Ice::Current&) {
	return conf->getReal(item);
}

/**
 * Get the Unit
 */
string	RCConfigI::getUnit(const string& item, const Ice::Current&) {
	return conf->getUnit(item);
}

/**
 * Get the Name
 */
string	RCConfigI::getName(const string& item, const Ice::Current&) {
	return conf->getName(item);
}

/**
 * Get the Index
 */
string	RCConfigI::getIndex(const string& item, const Ice::Current&) {
	return conf->getIndex(item);
}

/**
 *
 */
string	RCConfigI::getRange(const string& item, const Ice::Current&) {
	return conf->getRange(item);
}

/**
 * Get the Pre comment
 */
string	RCConfigI::getPreComment(const string& item, const Ice::Current&) {
	return conf->getPreComment(item);
}

/**
 * Get the Post comment
 */
string	RCConfigI::getPostComment(const string& item, const Ice::Current&) {
	return conf->getPostComment(item);
}

/**
 * Get the current Path from the Root level
 */
string	RCConfigI::getCurrentPath(const Ice::Current&) {
	return conf->getCurrentPath();
}

/**
 * Set an Int value
 */
void 	RCConfigI::setInt(const string& item, const string& value, const Ice::Current&) {
	conf->setInt(item, value);
}

/**
 * Set a String value
 */
void 	RCConfigI::setStr(const string& item, const string& value, const Ice::Current&) {
	conf->setStr(item, value);
}

/**
 * Set a Bool value
 */
void 	RCConfigI::setBool(const string& item, const string& value, const Ice::Current&) {
	conf->setBool(item, value);
}

/**
 * Set a Real value
 */
void 	RCConfigI::setReal(const string& item, const string& value, const Ice::Current&) {
	conf->setReal(item, value);
}

/**
 * Set the Unit
 */
void 	RCConfigI::setUnit(const string& item, const string& value, const Ice::Current&) {
	conf->setUnit(item, value);
}

/**
 * Set the Unit
 */
void 	RCConfigI::setName(const string& item, const string& value, const Ice::Current&) {
	conf->setName(item, value);
}

/**
 * Set of the Name
 */
void 	RCConfigI::setIndex(const string& item, const string& value, const Ice::Current&) {
	conf->setIndex(item, value);
}

/**
 * Set of the Index.
 */
void 	RCConfigI::setRange(const string& item, const string& value, const Ice::Current&) {
	conf->setRange(item, value);
}

/**
 * Set the Pre comment
 */
void 	RCConfigI::setPreComment(const string& item, const string& value, const Ice::Current&) {
	conf->setPreComment(item, value);
}

/**
 * Set the Post comment
 */
void 	RCConfigI::setPostComment(const string& item, const string& value, const Ice::Current&) {
	conf->setPostComment(item, value);
}

/**
 * Set path to "item" from root level
 */
void 	RCConfigI::setCurrentPath(const string& item, const Ice::Current&) {
	conf->setCurrentPath(item);
}


}
