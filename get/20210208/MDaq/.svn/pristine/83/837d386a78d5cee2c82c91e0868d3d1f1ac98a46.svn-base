/* =====================================================================================================================
 * CmdLineArgs.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::CmdLineArgs
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq data acquisition framework @ CEA-Irfu.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_utl_CmdLineArgs_h_INCLUDED
#define mdaq_utl_CmdLineArgs_h_INCLUDED

#include <vector>
#include <string>

namespace mdaq
{
namespace utl
{

/**
 * This is a helper class for going from an (argc,argv) paradigm to a string vector and vice-versa.
 */
class CmdLineArgs : public std::vector<std::string>
{
public:
	/**
	 * Builds a CmdLineArgs object from a given (argc, argv) pair.
	 * @param argc The argument count.
	 * @param argv The arguments c-vector
	 * @return The newly built CmdLineArg.
	 */
	static CmdLineArgs build(int argc, char** argv);

	/**
	 * Builds a CmdLineArgs object from a given function name a the string of its parameter arguments.
	 * @param funcName The c-string representing the calling function name.
	 * @param funcArgs The c-string containing all the command line arguments of the function call.
	 * @return The newly built CmdLineArg.
	 */
	static CmdLineArgs build(const char* funcName, const char* funcArgs);

public:
	/**
	 * Constructs a CmdLineArgs object. Reserves a minimum capacity of 128 for the argv vector.
	 * @param capacity Initial capacity of the argv vector.
	 */
	CmdLineArgs(size_t capacity = 128u) : argv_(0) { reserve(capacity); }

	CmdLineArgs(const std::vector< std::string > & args);

	/**
	 * Deletes allocated resources
	 */
	virtual ~CmdLineArgs() { destroy_argv(); }

	/**
	 * Appends a new argument at the end of the argument list.
	 * @param arg The argument to be appended.
	 * @return  The CmdLineArgs object itself so as to allow chained method calls.
	 */
	CmdLineArgs& addArg(const std::string& arg) { push_back(arg); return *this; }

	/**
	 * Erase all arguments.
	 * @return The CmdLineArgs object itself so as to allow chained method calls.
	 */
	CmdLineArgs& reset() { clear(); return *this; }

	/**
	 * Erases the i-th argument. If i is beyond the size of the argument list, nothing happens.
	 * @param i The index of the element to be erased.
	 * @return The CmdLineArgs object itself so as to allow chained method calls.
	 */
	CmdLineArgs& erase(size_type i);

	/**
	 * Returns an argv object corresponding to the argument list (first argument is a function name). The CmdLineArgs
	 * object keeps the ownership of the pointer.
	 * @return An argv pointer representing the argument list (including function name).
	 */
	char** argv();

	/**
	 * Returns the arguments count.
	 */
	int argc() const { return (int) size(); }

	/**
	 * Returns a reference to the i-th argument string.
	 * @param i The index of the argument (index 0 being the function name).
	 * @return The C++ string reference to the i-th argument.
	 */
	std::string& operator[] (size_type i) { return at(i); }

	/**
	 * Same as former operator, only for const access.
	 */
	const std::string& operator[] (size_type i) const { return at(i); }

private:
	void destroy_argv();

private:
	char** argv_; ///< This holds the argv. Its elements are re-computed at every call of argv().

};

} // namespace utl
} // namespace mdaq

#endif /* mdaq_utl_CmdLineArgs_h_INCLUDED */
