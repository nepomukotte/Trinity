/* =============================================================================
 * SymbolResolver.h
 * -----------------------------------------------------------------------------
 * class utl::stack::WinSymbolResolver
 * class utl::stack::DlAddrSymbolResolver
 * class utl::stack::VxSymbolResolver
 * Created on: 29 août 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef utl_stack_SymbolResolver_h_INCLUDED
#define utl_stack_SymbolResolver_h_INCLUDED

#if defined(_WIN32)
#   include <windows.h>
#   include <dbghelp.h>
#elif defined(linux)
#   include <dlfcn.h>
#elif defined(VXWORKS)
#   include <symLib.h>
typedef int intptr_t;
#endif

#ifndef DEMANGLER_BUFFER_SIZE
#define DEMANGLER_BUFFER_SIZE 256
#endif

namespace utl
{
namespace stack
{
/**
 * Default policy class to use when resolving the symbol name from a
 * stack-trace return address is either impossible or unneeded.
 */
template<class Base>
class NoSymbolResolver : public Base
{
};

#if defined(_WIN32)
/**
 * Policy class that resolves the symbol name from a stack-trace return address
 * using Windows API.
 */
template<class Base>
class WinSymbolResolver : public Base
{
public:
	WinSymbolResolver();

	const char* symbol();
	const char* mangledSymbol();
	intptr_t displacement();

private:
	const char* resolve();

	char buffer[sizeof(SYMBOL_INFO)+DEMANGLER_BUFFER_SIZE];
	DWORD64 disp;
};
#endif

#if defined(linux)

/**
 * Policy class that resolves the symbol name from a stack-trace return address
 * using dladdr unix function.
 */
template<class Base>
class DlAddrSymbolResolver : public Base
{
public:
	DlAddrSymbolResolver();
	~DlAddrSymbolResolver();

	const char* symbol();
	const char* mangledSymbol();
	const char* module();
	intptr_t displacement();

private:
	void resolve();
	void demangle();

	Dl_info dlinfo;
	char* demangledSymbol;
};
#endif

#if defined(VXWORKS)

/**
 * Policy class that resolves the symbol name from a stack-trace return address
 * using VxWorks API.
 */
template<class Base>
class VxSymbolResolver : public Base
{
public:
	VxSymbolResolver();
	~VxSymbolResolver();

	const char* symbol();
	const char* mangledSymbol();
	intptr_t displacement();

private:
	void resolve();
	void demangle();

	char* symbolName;
	int symbolValue;
	SYM_TYPE symbolType;
	char demangledSymbol[DEMANGLER_BUFFER_SIZE];
};

#endif

} /* namespace stack */
} /* namespace utl */

#endif /* utl_stack_SymbolResolver_h_INCLUDED */
