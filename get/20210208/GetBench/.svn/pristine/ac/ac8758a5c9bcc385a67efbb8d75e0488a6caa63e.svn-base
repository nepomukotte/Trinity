/**
 * @file Node.h
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

#ifndef get_rc_Node_h_INCLUDED
#define get_rc_Node_h_INCLUDED

#include "CCfg/CConfig.h"
#include "get/rc/DataLinkSet.h"
#include "get/rc/DataSenderSet.h"
#include "get/rc/NodeId.h"
#include "mdaq/rc/HardwareNode.h"
#include "mdaq/hw/Control.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/exception_ptr.hpp>

//_________________________________________________________________________________________________
namespace get
{
namespace rc
{

class BackEnd;

/**
 * This abstract class encapsulates the connection to an network endpoint containing a hardware object such as a CoBo
 * or a Mutant. It is expected to be specialized by specific objects, e.g. CoBoNode or MutantNode.
 * Every Node instance has a unique NodeId (e.g. CoBo[3]) and is initialized through a CConfig object.
 * The main function of this class is to respond to the global state machine transitions; specialized classes such as
 * CoBoNode must essentially implement the pure virtual methods that respond to state transitions (impl_start(),
 * impl_stop(), etc.).
 */
class Node
{
public:
	class List;
	class ListList;

public:
	Node(const NodeId& nodeId, get::rc::BackEnd* backend=0, mdaq::rc::HardwareNode* _hwNode=0);
	virtual ~Node() {}

	static CCfg::CConfig nodeCfg(CCfg::CConfig globalCfg, const NodeId& nodeId)
	{
		return globalCfg("Node", nodeId.typeId)("Instance", nodeId.instanceId);
	}

	CCfg::CConfig nodeCfg(CCfg::CConfig globalCfg) const { return nodeCfg(globalCfg, id()); }

	const NodeId& id() const { return nodeId; }

	virtual void setEnabled(bool enabled) { enabled_ = enabled; }
	virtual bool isEnabled() const { return enabled_; };

	virtual DataSenderSet dataSenders() const = 0;
	virtual DataLinkSet dataLinks() const;
	virtual void setDataLinks(const DataLinkSet & links);
	virtual void daqConnect();
	virtual void daqDisconnect();

	const mdaq::rc::HardwareNode & hwNode() const { return *hwNode_; };
	      mdaq::rc::HardwareNode & hwNode()       { return *hwNode_; };
public:
	void call_prepare   (CCfg::CConfig nodeCfg) { impl_prepare(nodeCfg)  ; }
	void call_unprepare (CCfg::CConfig        ) { impl_unprepare()       ; }
	void call_configure (CCfg::CConfig nodeCfg) { impl_configure(nodeCfg); }
	void call_start     (CCfg::CConfig        ) { impl_start()           ; }
	void call_stop      (CCfg::CConfig        ) { impl_stop()            ; }

protected:
	virtual void impl_prepare   (CCfg::CConfig nodeCfg) = 0;
	virtual void impl_unprepare (                     ) = 0;
	virtual void impl_configure (CCfg::CConfig nodeCfg) = 0;
	virtual void impl_start     (                     ) = 0;
	virtual void impl_stop      (                     ) = 0;

protected:
	NodeId nodeId;
	mdaq::hw::DevicePrx device(const std::string & name) { return hwNode().findDevice(name); }
	get::rc::BackEnd* backend;
	mdaq::rc::HardwareNode* hwNode_;
	DataLinkSet dataLinks_;
	bool enabled_;
};

#define NODE_DEBUG() LOG_DEBUG() << id().toString() << ": "
#define NODE_INFO()  LOG_INFO()  << id().toString() << ": "
#define NODE_WARN()  LOG_WARN() << id().toString() << ": "
#define NODE_ERROR() LOG_ERROR() << id().toString() << ": "

/**
 * This class embodies a list of Node instances of a certain type (e.g. a list of CoBo nodes). The purpose of this
 * container is to implement operations that execute the same function on all contained nodes (see doLoop methods
 * below). These executions can be accomplished either in parallel or in sequence.
 */
class Node::List
{
public:
	friend class Node;
	typedef boost::ptr_map<std::string, Node> NodeMap;
	typedef std::auto_ptr<Node> (*NodeAllocator) (const std::string&, mdaq::rc::HardwareNode*, get::rc::BackEnd*);
	typedef std::map<std::string, NodeAllocator> NodeFactory;

	/**
	 * This template class is an iterator over all the elements of Node::List. The iteration is "Java-like" as it uses
	 * a hasNext() logic instead of a the begin()-end() style common in C++ containers.
	 * @tparam C Template parameter excpected to specify the actual specialized class of the contained Node (e.g.
	 *           CoBoNode)
	 */
	template<class C>
	class Iterator
	{
	public:
		Iterator(NodeMap& nodeMap) : nodeMap(nodeMap), it(nodeMap.begin()) {}
		bool hasNext() const { return it != nodeMap.end(); }
		C& node() { return dynamic_cast<C&>(*(it->second)); }
		C& operator *() { return node(); }
		C* operator->() { return dynamic_cast<C*>(it->second); }
		std::string key() const { return it->first; }
		Iterator& operator ++() { ++it; return *this; }

	private:
		NodeMap& nodeMap;
		NodeMap::iterator it;
	};

	/**
	 * Same as ietrator above but for const access to Nodes.
	 */
	template<class C>
	class ConstIterator
	{
	public:
		ConstIterator(const NodeMap& nodeMap) : nodeMap(nodeMap), it(nodeMap.begin()) {}
		bool hasNext() const { return it != nodeMap.end(); }
		const C& node() const { return dynamic_cast<const C&>(*(it->second)); }
		const C& operator *() const { return node(); }
		const C* operator->() const { return dynamic_cast<const C*>(it->second); }
		std::string key() const { return it->first; }
		ConstIterator& operator ++() { ++it; return *this; }

	private:
		const NodeMap& nodeMap;
		NodeMap::const_iterator it;
	};

	/**
	 * The values of this enum are used to specify a an operation on all nodes of the list must be executed in parallel
	 * or in sequence.
	 */
	enum ExecMode { InSequence, InParallel };

public:
	List(const std::string&  nodeType);
	size_t size() const { return nodeMap.size(); } ///< @return The number of nodes contained by this Node::List.

	template<class C> Iterator<C> iterate();
	template<class C> ConstIterator<C> iterate() const;
	template<class C> static std::auto_ptr<Node>
		nodeAllocate(const std::string& instanceId, mdaq::rc::HardwareNode* hwNode, get::rc::BackEnd* backend);

	static CCfg::CConfig nodeCfg(CCfg::CConfig globalCfg, const NodeId& nodeId);

	List& addNodeAllocator(const std::string& nodeType, NodeAllocator nodeAllocator)
	{
		nodeFactory[nodeType] = nodeAllocator;
		return *this;
	}

	List& insertNode(const std::string& instanceId, std::auto_ptr<Node> newNode)
	{
		nodeMap.insert(instanceId, newNode);
		return *this;
	}

	List& insertNewNode (
		const std::string& typeId,
		const std::string& instanceId,
		mdaq::rc::HardwareNode* hwNode=0,
		get::rc::BackEnd* backend=0
	)
	{
		return insertNode(instanceId, nodeFactory[typeId](instanceId, hwNode, backend));
	}

	      Node& getNode(const std::string& instanceId)       { try { return nodeMap.at(instanceId); } catch (const boost::bad_ptr_container_operation & ) { throw mdaq::utl::CmdException(std::string("Could not find node instance named '") + instanceId + "'.");}}
	const Node& getNode(const std::string& instanceId) const {  try { return nodeMap.at(instanceId); } catch (const boost::bad_ptr_container_operation & ) { throw mdaq::utl::CmdException(std::string("Could not find node instance named '") + instanceId + "'.");}}
	template <class C>       C& getNode(const std::string& instanceId)       { return dynamic_cast<C&>(getNode(instanceId)); }
	template <class C> const C& getNode(const std::string& instanceId) const { return dynamic_cast<C&>(getNode(instanceId)); }

	      Node& operator[] (const std::string& instanceId)       { return getNode(instanceId); }
	const Node& operator[] (const std::string& instanceId) const { return getNode(instanceId); }

	void removeNodes() { nodeMap.clear(); }

	void doLoop(ExecMode execMode, void (Node::*method) (CCfg::CConfig), CCfg::CConfig& globalCCfg);
	void doLoop(void (Node::*method) (CCfg::CConfig), CCfg::CConfig& globalCCfg)
	{
		doLoop(globalCCfg("ECC")("parallelProcessing").getBoolValue(true)?InParallel:InSequence, method, globalCCfg);
	}
	void doLoop_prepare   (CCfg::CConfig globalCCfg) { doLoop(&Node::call_prepare  , globalCCfg); }
	void doLoop_unprepare (CCfg::CConfig globalCCfg) { doLoop(InSequence, &Node::call_unprepare, globalCCfg); }
	void doLoop_configure (CCfg::CConfig globalCCfg) { doLoop(&Node::call_configure, globalCCfg); }
	void doLoop_start     (CCfg::CConfig globalCCfg) { doLoop(&Node::call_start    , globalCCfg); }
	void doLoop_stop      (CCfg::CConfig globalCCfg) { doLoop(InSequence, &Node::call_stop     , globalCCfg); }

private:
	static void call_method_and_catch_exceptions(void (Node::*method) (CCfg::CConfig), Node* node, CCfg::CConfig& globalCfg, boost::exception_ptr & error);

public:
	const std::string nodeType;

private:
	NodeMap nodeMap;
	NodeFactory nodeFactory;
};

/**
 * Template function returning an iterator over the specialized Node classes of this Node::List.
 * @return The iterator instantiated with the actual class of the specialized Node.
 */
template<class C>
Node::List::Iterator<C>
Node::List::iterate()
{
	return Iterator<C>(nodeMap);
}

/**
 * Same as previous template function for the const iterator.
 */
template<class C>
Node::List::ConstIterator<C>
Node::List::iterate() const
{
	return ConstIterator<C>(nodeMap);
}

/**
 * Allocates a new specialized concrete Node.
 * @tparam C The class of the concrete Node to be instantiated.
 * @param instanceId The ID of this new node instance.
 * @param hwNode The hardware node associated with this node (@see HardWareNode).
 * @param backend A pointer to the Backend object containing this Node::List.
 * @return The newly allocated concrete Node.
 */
template<class C>
std::auto_ptr<Node> Node::List::nodeAllocate(const std::string& instanceId, mdaq::rc::HardwareNode* hwNode, get::rc::BackEnd* backend)
{
	return std::auto_ptr<Node>(new C(instanceId, hwNode, backend));
}

class Node::ListList
{
public:
	typedef boost::ptr_map<std::string, List> ListMap;

public:
	ListList(BackEnd& backend_) : backend_(backend_) {}

	List& addList(std::auto_ptr<List> newListPtr)
	{
		const std::string nodeType(newListPtr->nodeType); // Make a copy of this because auto_ptr will be deallocated.
		return *(listMap.insert(nodeType, newListPtr).first->second);
	}

	template <class C>
	List& addNewList()
	{
		List& newList = addList(std::auto_ptr<List>(new Node::List(C::NODE_TYPE)));
		newList.addNodeAllocator(C::NODE_TYPE, &List::nodeAllocate<C>);
		return newList;
	}

	// @TODO: replace addNewList by putting instantiation in getList (instantiates if not already there)
	      Node::List& getList(const std::string& nodeType)       { return listMap.at(nodeType); }
	const Node::List& getList(const std::string& nodeType) const { return listMap.at(nodeType); }

	template <class C>
	Node::List& getList() { return getList(C::NODE_TYPE); }

	template <class C>
	const Node::List& getList() const { return getList(C::NODE_TYPE); }

	void addNewNode(const std::string& nodeType, const std::string& instanceId, ::mdaq::rc::HardwareNode& hwNode)
	{
		listMap.at(nodeType).insertNewNode(nodeType, instanceId, &hwNode, &backend_);
	}

	template <class C>
	void addNewNode(const std::string& instanceId, ::mdaq::rc::HardwareNode& hwNode)
	{
		addNewNode(C::NODE_TYPE, instanceId, hwNode);
	}

	Node& getNode(const std::string& nodeType, const std::string& instanceId)
	{
		return getList(nodeType).getNode(instanceId);
	}

	const Node& getNode(const std::string& nodeType, const std::string& instanceId) const
	{
		return getList(nodeType).getNode(instanceId);
	}

	template <class C>       C& getNode(const std::string& instanceId)       { return dynamic_cast<      C&>(getNode(C::NODE_TYPE, instanceId)); }
	template <class C> const C& getNode(const std::string& instanceId) const { return dynamic_cast<const C&>(getNode(C::NODE_TYPE, instanceId)); }

	void removeNodes();

	template <class C> Node::List::     Iterator<C> iterate()       { return getList<C>().template iterate<C>(); }
	template <class C> Node::List::ConstIterator<C> iterate() const { return getList<C>().template iterate<C>(); }

protected:
	std::auto_ptr<mdaq::hw::RemoteRegisterProcessor> regProcPtr;

private:
	mutable ListMap listMap;
	BackEnd& backend_;
};

//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_Node_h_INCLUDED */
