/* =============================================================================
 * SmartRegister.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::SmartRegister
 * Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of
 * the Mordicus Real-Time Software Development Framework.
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

#ifndef mdaq_hw_SmartRegister_hpp_INCLUDED
#define mdaq_hw_SmartRegister_hpp_INCLUDED

#include "AbstractRegister.hpp"
#include "NestedStoragePolicy.hpp"
#include "ObservableValue.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace mdaq
{
namespace hw
{
/**
 * This is a policy-based class representing a hardware register.
 * @tparam T The datatype of the register. It can be either primitive types like
 * uint_32 and uint_8, or classes (to implement non-standard register sizes)
 * having assignment and arithmetic operators defined.
 * @tparam S The storage policy: it implements access to the register and
 * contains any necessary data to do this. The policy class must define the
 * following 2 methods that access the register.
 * <ul>
 *    <li>T read() const
 *    <li>void write(T value)
 * </ul>
 * It can also define the following 2 typedefs, if the register is addressable
 * from the CPU:
 * <ul>
 *    <li>AddressType
 * </ul>
 * @tparam O The observation policy: it defines how modifications to the
 * register can be observed by an outside actor.
 * @tparam B The base class policy: it can be used to enable polymorphism or
 * to introduce additionnal methods or/and data.
 */
template <typename T,
	template <class Type> class S = NestedStoragePolicy,
	template <class Type> class O = NonObservableValue,
	class B = DefaultBaseRegister>
class SmartRegister : public B, public S<T>, public O<T>
{
public:
	typedef T WordType;			///< The word type of the register
	typedef S<T> StoragePolicy;	///< How the register is stored and accessed.
	typedef O<T> ObservationPolicy;///< How the register can be observed.
	typedef B    BaseRegister;		///< any other base class (may be used to introduce polymorphism...).

	/// Constructor.
	SmartRegister() {}
	SmartRegister(T value) : StoragePolicy(value) {}
	SmartRegister(typename StoragePolicy::AddressType address) : StoragePolicy(address) {}
	SmartRegister(typename StoragePolicy::DeviceType* device, typename StoragePolicy::AddressType address) : StoragePolicy(device, address) {}

	/// Sets the register to specified value.
	/// @param val value written to registers using access policy
	void setValue(WordType val) { observableWrite<SmartRegister>(*this, val); }

	/// The value of the register.
	/// @return Register's value as read by the access policy
	WordType value() const { return StoragePolicy::read(); }

	/// Cast operator to the register's word type.
	operator WordType() { return value(); }

	/// @name All register operators called using the access policy
	//@{
	SmartRegister& operator   = (WordType val) { setValue(val); return *this;   }
	SmartRegister& operator  &= (WordType val) { return *this = value() &  val; }
	SmartRegister& operator  |= (WordType val) { return *this = value() |  val; }
	SmartRegister& operator  ^= (WordType val) { return *this = value() ^  val; }
	SmartRegister& operator <<= (int val)      { return *this = value() << val; }
	SmartRegister& operator >>= (int val)      { return *this = value() >> val; }
	SmartRegister& operator  += (WordType val) { return *this = value() +  val; }
	SmartRegister& operator  -= (WordType val) { return *this = value() -  val; }
	SmartRegister& operator  *= (WordType val) { return *this = value() *  val; }
	SmartRegister& operator  /= (WordType val) { return *this = value() /  val; }
	SmartRegister& operator  %= (WordType val) { return *this = value() %  val; }
	SmartRegister& operator  ++ ()		{ return *this += 1; }
	SmartRegister& operator  -- ()		{ return *this -= 1; }
	WordType       operator  ++ (int)	{ WordType w = *this; ++(*this); return w; }
	WordType       operator  -- (int)	{ WordType w = *this; --(*this); return w; }
	//@}

private:
	/// Writes to the register and notify observers.
	/// This overloaded function is disabled (using the SFINAE idiom) if the class is not using an observation policy.
	template<class R>
	void observableWrite(typename boost::disable_if<boost::is_base_of<NonObservableValue<T>, R>, R>::type& obj, T val)
	{
		T oldValue = obj.value();
		obj.write(val);
		obj.notifyObservers(oldValue, val);
	}

	/// Write to the register without notifying observers.
	/// This overloaded function is enabled (using the SFINAE idiom) if the class is not using an observation policy.
	template<class R>
	void observableWrite(typename boost::enable_if<boost::is_base_of<NonObservableValue<T>, R>, R>::type& obj, T val)
	{
		obj.write(val);
	}

};

} // namespace hw
} // namespace mdaq

#endif // mdaq_hw_SmartRegister_hpp_INCLUDED
