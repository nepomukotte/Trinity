/**
 * @file Node.cpp
 * @date 15 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "get/rc/Node.h"
#include <utl/Logging.h>

using mdaq::utl::CmdException;

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
Node::Node(const NodeId& nodeId, get::rc::BackEnd* backend, mdaq::rc::HardwareNode* node)
	: nodeId(nodeId), backend(backend), hwNode_(node), enabled_(true)
{
}
//_________________________________________________________________________________________________
DataLinkSet Node::dataLinks() const
{
	return dataLinks_;
}
//_________________________________________________________________________________________________
/**
 * Stores list of data routers related to this node's data senders.
 */
void Node::setDataLinks(const DataLinkSet & links)
{
	dataLinks_.clear();
	DataSenderSet senders = dataSenders();
	DataLinkSet::const_iterator linkIt;
	for (linkIt = links.begin(); linkIt != links.end(); ++linkIt)
	{
		if (find(senders.begin(), senders.end(), (*linkIt).sender) != senders.end())
			dataLinks_.insert(*linkIt);
	}
	LOG_DEBUG() << "Registered " << dataLinks_.size() << " data link(s) for node " << hwNode().hwServerCtrlPrx->name();
}
//_________________________________________________________________________________________________
void Node::daqConnect()
{
	DataSenderSet senders = dataSenders();
	DataSenderSet::const_iterator senderIter;
	for (senderIter = senders.begin(); senderIter != senders.end(); ++ senderIter)
	{
		// Find link for sender
		DataLinkSet::const_iterator linkIter;
		linkIter = std::find_if(dataLinks_.begin(), dataLinks_.end(), IsLinkSender(*senderIter));
		if (linkIter != dataLinks_.end())
		{
			try
			{
				LOG_INFO() << "Connecting data sender '" << (*linkIter).sender.id << "' to router '"
						<< (*linkIter).router.name << "' of type '" << (*linkIter).router.flowType
						<< "' at address " << (*linkIter).router.sockAddr;
				hwNode().connectRunProcessor((*linkIter).router.sockAddr.toString(),
						(*linkIter).router.flowType);
			}
			catch (const CmdException & e)
			{
				throw CmdException(std::string("Could not connect data sender '") + (*linkIter).sender.id
								+ "' to its data router: " + e.reason);
			}
			catch (const std::exception & e)
			{
				throw CmdException(std::string("Could not connect data sender '") + (*linkIter).sender.id
						+ "' to its data router: " + e.what());
			}
		}
		else
		{
			throw CmdException(std::string("Could not find a valid data router for data sender '") + (*senderIter).id + '\'');
		}
	}
}
//_________________________________________________________________________________________________
/**
 * Disconnects from data routers.
 */
void Node::daqDisconnect()
{
	try
	{
		hwNode().disconnectRunProcessor();
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << e.reason;
	}
	catch (const std::exception & e)
	{
		LOG_WARN() << e.what();
	}
}

void Node::List::call_method_and_catch_exceptions(void (Node::*method) (CCfg::CConfig), Node* node, CCfg::CConfig& nodeCfg, boost::exception_ptr & error)
{
	try
	{
		try
		{
			(node->*method)(nodeCfg);
		}
		catch (const mdaq::utl::CmdException & e)
		{
			throw boost::enable_current_exception(e);
		}
		error = boost::exception_ptr();
	}
	catch (...)
	{
		error = boost::current_exception();
	}
}

void Node::List::doLoop (
	ExecMode execMode,
	void (Node::*method) (CCfg::CConfig),
	CCfg::CConfig& globalCfg
)
{
	bool hasConfig = globalCfg.getAttributes().size() > 0;
	switch(execMode)
	{
	case InSequence:
		// Loop over contained objects
		for (NodeMap::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it)
		{
			(it->second->* method)(hasConfig?nodeCfg(globalCfg, it->second->id()):globalCfg);
		}
		break;

	case InParallel:
		// Prepare thread group for all contained objects
		boost::thread_group objectThreads;

		// Launch one thread per contained objects
		typedef std::vector< boost::exception_ptr > ErrorList;
		ErrorList errors(nodeMap.size());
		ErrorList::iterator errorIt = errors.begin();
		for (NodeMap::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it, ++errorIt)
		{
			objectThreads.create_thread(
					boost::bind(
							&Node::List::call_method_and_catch_exceptions,
							method,
							it->second,
							hasConfig?nodeCfg(globalCfg, it->second->id()):globalCfg,
							boost::ref(*errorIt)));
		}

		// Wait for all parallel threads to end
		objectThreads.join_all();

		// Handle exceptions
		std::ostringstream errorMsg;
		errorIt = errors.begin();
		for (NodeMap::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it, ++errorIt)
		{
			if (*errorIt)
			{
				try
				{
					boost::rethrow_exception(*errorIt);
				}
				catch (CmdException & e)
				{
					errorMsg << it->second->id().toString() << ": " << e.reason << " | ";
				}
			}
		}
		if (not errorMsg.str().empty())
		{
			throw mdaq::utl::CmdException(errorMsg.str());
		}
		break;
	}
}

void Node::ListList::removeNodes()
{
	for (ListMap::iterator itList = listMap.begin(); itList != listMap.end(); ++itList)
	{
		itList->second->removeNodes();
	}
}

/**
 * This is the constructor of the Node::List, which is expected to contain instances of a single type of specialized
 * Node.
 * @param nodeType The string representation of the type of nodes contained in this List (e.g. "CoBo").
 */
Node::List::List(const std::string&  nodeType) : nodeType(nodeType)
{
	LOG_DEBUG() << "--- Construct Node::List of type: " << nodeType;
}

/**
 * This static function extracts from a global CConfig the sub-config relating to a specific Node type and instance
 * @param globalCfg The global CConfig object.
 * @param nodeId The ID os the specialized node specifying its type and instance ids.
 * @return The sub-config of the specified node.
 */
CCfg::CConfig
Node::List::nodeCfg(CCfg::CConfig globalCfg, const NodeId& nodeId)
{
	return globalCfg("Node", nodeId.typeId)("Instance", nodeId.instanceId);
}


} /* namespace rc */
} /* namespace get */
