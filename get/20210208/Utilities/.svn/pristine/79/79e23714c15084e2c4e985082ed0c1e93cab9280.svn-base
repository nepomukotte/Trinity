/* =============================================================================
 * Singleton.ipp
 * -----------------------------------------------------------------------------
 * class utl::Singleton
 * Created on: 3 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities project.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifndef utl_Singleton_ipp_INCLUDED
#define utl_Singleton_ipp_INCLUDED

#include "Singleton.hpp"
#include <cstdlib>

namespace utl
{
/**
 * Indicates whether the instance has been destroyed as part of program shut down.
 * This is necessary because static "by value" storage gives no clue about that.
 */
template<class T>
bool Singleton<T>::destroyed = false;

/**
 * Points to the current instance of the singleton.
 */
template<class T>
Singleton<T>* Singleton<T>::instance_ = 0;

/**
 * Constructor.
 * Initializes the T element.
 */
template<class T>
Singleton<T>::Singleton()
{
}

/**
 * Destructor.
 */
template<class T>
Singleton<T>::~Singleton()
{
	destroyed = true;
	instance_ = 0;
}

/**
 * Resurrect the singleton when instance() is called during global destructors
 * execution by another singleton that outlives this one.
 */
template<class T>
void Singleton<T>::resurrect()
{
	new(instance_) Singleton<T>;
	std::atexit(destroy);
	destroyed = false;
}

/**
 * Destroys the singleton (callback for std::atexit).
 */
template<class T>
void Singleton<T>::destroy()
{
	instance_->~Singleton<T>();
}

/**
 * Return the singleton instance.
 */
template<class T>
T& Singleton<T>::instance()
{
	if(!instance_)
	{
		static Singleton<T> theInstance;
		if(destroyed)
		{
			resurrect();
		}
		instance_ = &theInstance;
	}
	return instance_->element;
}

} /* namespace utl */

#endif /* utl_Singleton_ipp_INCLUDED */
