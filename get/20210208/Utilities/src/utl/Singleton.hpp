/* =============================================================================
 * Singleton.hpp
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

#ifndef utl_Singleton_hpp_INCLUDED
#define utl_Singleton_hpp_INCLUDED

namespace utl
{
/**
 * This is a generic singleton class handling object lifetime (dead reference
 * problem) correctly.
 * It inspires from "Modern C++ Design" by Andrei Alexandrescu. Although
 * currently no policy is exposed, it may be improved in the future to support
 * this kind of customization.
 */
template<class T>
class Singleton
{
public:
	typedef T elem_type;

	static T& instance();

private:
	Singleton();
	Singleton(const Singleton&);
	~Singleton();
	Singleton& operator=(const Singleton&);

	T element;

	static void destroy();
	static void resurrect();

	static Singleton<T>* instance_;
	static bool destroyed;
};

/**
 * This class is made to wrap a class that may be instanciated as a global
 * several times, while ensuring its constructor is only executed once.
 * It is mainly useful to perform some initialization as a side-effect of
 * construction of globals, before main().
 */
template<class InitClass>
struct Once
{
	Once()
	{
		static InitClass doInit;
	}
};

} /* namespace utl */

#endif /* utl_Singleton_hpp_INCLUDED */
