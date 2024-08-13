/**
 * @file NodeMapLoop.hpp
 * @date Aug 25, 2016
 * @author anvar
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Shebli.Anvar@cea.fr
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

#ifndef get_utl_NodeMapLoop_hpp_INCLUDED
#define get_utl_NodeMapLoop_hpp_INCLUDED

#include "CCfg/CConfig.h"
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <map>

namespace get {
namespace utl {

enum ExecMode { InSequence, InParallel };

template <typename K, typename T, ExecMode E>
struct NodeMapLoop;

template <typename K, typename T>
struct NodeMapLoop<K, T, InSequence>
{
	typedef K KeyType;
	typedef T ObjectType;
	typedef typename boost::ptr_map<KeyType, ObjectType> ContainerType;

	static void call(ContainerType& container, void (ObjectType::*method)(CCfg::CConfig ccfg))
	{
		// Loop over contained objects
		for (typename ContainerType::iterator it = container.begin(); it != container.end(); ++it)
		{
			it->second.*method(ccfg("Node", id()));
		}
	}
};

template <typename K, typename T>
struct NodeMapLoop<K, T, InParallel>
{
	typedef K KeyType;
	typedef T ObjectType;
	typedef typename boost::ptr_map<KeyType, ObjectType> ContainerType;

	static void call(ContainerType& container, void (ObjectType::*method)())
	{
		// Prepare thread group for all contained objects
		boost::thread_group objectThreads;

		// Launch one thread per contained objects
		for (typename ContainerType::iterator it = container.begin(); it != container.end(); ++it)
		{
			objectThreads.create_thread(boost::bind(method, &(it->second)));
		}

		// Wait for all parallel threads to end
		objectThreads.join_all();
	}
};

} // namespace utl
} // namespace get



#endif /* get_utl_NodeMapLoop_hpp_INCLUDED */
