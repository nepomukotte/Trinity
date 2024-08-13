/*
 * $Id: VersaPointer.hpp 1780 2017-02-10 8:53:16Z anvar $
 * @file VersaPointer.hpp
 * @created Feb 10, 2017
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */


#ifndef mdaq_utl_VersaPointer_hpp_INCLUDED
#define mdaq_utl_VersaPointer_hpp_INCLUDED

#include "utl/Logging.h"
#include <boost/noncopyable.hpp>
#include <memory>

namespace mdaq
{
namespace utl
{

template <class C>
class VersaPointer : boost::noncopyable
{
public:
	typedef C element_type;

public:
	VersaPointer(element_type* newElement);
	VersaPointer(std::auto_ptr<element_type> newElement);
	VersaPointer(element_type& extElement);
	~VersaPointer();
	void reset(element_type* newElement);
	void reset(std::auto_ptr<element_type> newElement);
	void reset(element_type& newElement);
	bool owned() { return ownedElement.get() != 0; }
	const element_type* get         () const { return  p_element; }
	      element_type* get         ()       { return  p_element; }
	const element_type& element     () const { return *p_element; }
	      element_type& element     ()       { return *p_element; }
	const element_type* operator -> () const { return  p_element; }
	      element_type* operator -> ()       { return  p_element; }
	const element_type& operator *  () const { return *p_element; }
	      element_type& operator *  ()       { return *p_element; }

private:
	element_type* p_element;
	std::auto_ptr<C> ownedElement;
};

template <class C>
VersaPointer<C>::VersaPointer(element_type* newElement) : p_element(0)
{
	reset(newElement);
}

template <class C>
VersaPointer<C>::VersaPointer(std::auto_ptr<element_type> newElement) : p_element(0)
{
	reset(newElement.release());
}

template <class C>
VersaPointer<C>::VersaPointer(element_type& extElement) : p_element(0)
{
	reset(extElement);
}

template <class C>
void VersaPointer<C>::reset(element_type* newElement)
{
	p_element = newElement;
	ownedElement.reset(p_element);
}

template <class C>
void VersaPointer<C>::reset(std::auto_ptr<element_type> newElement)
{
	p_element = newElement.release();
	ownedElement.reset(p_element);
}

template <class C>
void VersaPointer<C>::reset(element_type& extElement)
{
	p_element = &extElement;
	ownedElement.reset();
}

template <class C>
VersaPointer<C>::~VersaPointer()
{
	ownedElement.reset();
}

} /* namespace utl */
} /* namespace mdaq */

#endif /* mdaq_utl_VersaPointer_hpp_INCLUDED */
