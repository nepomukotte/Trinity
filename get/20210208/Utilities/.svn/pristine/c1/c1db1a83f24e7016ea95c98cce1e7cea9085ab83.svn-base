/* =============================================================================
 * SymbolResolver.ipp
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

#include "SymbolResolver.h"

#if defined(linux)
#   include <cxxabi.h>
#elif defined(VXWORKS)
#   include <vxWorks.h>
#   include <sysSymTbl.h>
#   include <demangler.h>
#endif

namespace utl
{
namespace stack
{

#if defined(_WIN32)

/**
 * Default constructor.
 */
template<class Base>
WinSymbolResolver<Base>::WinSymbolResolver()
	: disp(0xFFFFFFFFFFFFFFFF)
{
}

/**
 * Performs the symbol resolution for the address returned by Base::address().
 * This method is designed to do its job only once per object, using
 * lazy-initialization.
 */
template<class Base>
const char* WinSymbolResolver<Base>::resolve()
{
	SYMBOL_INFO* symInfo = reinterpret_cast<SYMBOL_INFO*>(buffer);
	if(disp == 0xFFFFFFFFFFFFFFFF)
	{
		symInfo->SizeOfStruct = sizeof(SYMBOL_INFO);
		symInfo->MaxNameLen = sizeof(buffer) - sizeof(SYMBOL_INFO) + 1;
		disp = 0;
		SymFromAddr(GetCurrentProcess(), reinterpret_cast<DWORD64>(this->address()), &disp, symInfo);
	}
	return reinterpret_cast<const char*>(&symInfo->Name);
}

/**
 * Returns the demangled symbol name.
 */
template<class Base>
const char* WinSymbolResolver<Base>::symbol()
{
	return resolve();
}

/**
 * Returns the raw (possibly mangled) symbol name.
 */
template<class Base>
const char* WinSymbolResolver<Base>::mangledSymbol()
{
	return resolve();
}

/**
 * Return the offset in bytes between the function beginning (in machine code)
 * and the stack-frame return address.
 */
template<class Base>
intptr_t WinSymbolResolver<Base>::displacement()
{
	resolve();
	return static_cast<intptr_t>(disp);
}
#endif // defined(_WIN32)



#if defined(linux)

/**
 * Default constructor.
 */
template<class Base>
DlAddrSymbolResolver<Base>::DlAddrSymbolResolver()
	: demangledSymbol(0)
{
	dlinfo.dli_fbase = 0;
}

/**
 * Destructor: free allocated memory for symbol demangling.
 */
template<class Base>
DlAddrSymbolResolver<Base>::~DlAddrSymbolResolver()
{
	free(demangledSymbol);
}

/**
 * Performs the symbol resolution for the address returned by Base::address().
 * This method is designed to do its job only once per object, using
 * lazy-initialization.
 */
template<class Base>
void DlAddrSymbolResolver<Base>::resolve()
{
	if(dlinfo.dli_fbase == 0)
	{
		dladdr(this->address(), &dlinfo);
		if(dlinfo.dli_sname == 0)
		{
			dlinfo.dli_sname = "";
			dlinfo.dli_saddr = this->address();
		}
		if(dlinfo.dli_fname == 0)
		{
			dlinfo.dli_fname = "";
		}
	}
}

/**
 * Performs demangling of the previously resolved symbol.
 */
template<class Base>
void DlAddrSymbolResolver<Base>::demangle()
{
	if(demangledSymbol == 0)
	{
		size_t length = 0;
		int status;
		char* ret = abi::__cxa_demangle(dlinfo.dli_sname, demangledSymbol, &length, &status);
		if(status == 0)
		{
			demangledSymbol = ret;
		}
		else
		{
			free(ret);
		}
	}
}

/**
 * Returns the demangled symbol name.
 */
template<class Base>
const char* DlAddrSymbolResolver<Base>::symbol()
{
	resolve();
	demangle();
	return demangledSymbol != 0 ? demangledSymbol : dlinfo.dli_sname;
}

/**
 * Returns the raw (possibly mangled) symbol name.
 */
template<class Base>
const char* DlAddrSymbolResolver<Base>::mangledSymbol()
{
	resolve();
	return dlinfo.dli_sname;
}

/**
 * Return the full path of the module where the symbol is defined.
 */
template<class Base>
const char* DlAddrSymbolResolver<Base>::module()
{
	resolve();
	return dlinfo.dli_fname;
}

/**
 * Return the offset in bytes between the function beginning (in machine code)
 * and the stack-frame return address.
 */
template<class Base>
intptr_t DlAddrSymbolResolver<Base>::displacement()
{
	resolve();
	return reinterpret_cast<intptr_t>(this->address())
		   - reinterpret_cast<intptr_t>(dlinfo.dli_saddr);
}

#endif // defined(linux)




#if defined(VXWORKS)

/**
 * Default constructor.
 */
template<class Base>
VxSymbolResolver<Base>::VxSymbolResolver()
	: symbolName(0), symbolValue(-1)
{
	demangledSymbol[0] = '\0';
}

/**
 * Destructor: free allocated memory for symbol resolution.
 */
template<class Base>
VxSymbolResolver<Base>::~VxSymbolResolver()
{
	free(symbolName);
}

/**
 * Performs the symbol resolution for the address returned by Base::address().
 * This method is designed to do its job only once per object, using
 * lazy-initialization.
 */
template<class Base>
void VxSymbolResolver<Base>::resolve()
{
	if(symbolValue == -1)
	{
		if(symByValueFind(sysSymTbl,
			reinterpret_cast<int>(this->address()),
			&symbolName, &symbolValue, &symbolType)  != OK)
		{
			symbolName = 0;
		}
	}
}

/**
 * Performs demangling of the previously resolved symbol.
 */
template<class Base>
void VxSymbolResolver<Base>::demangle()
{
	if(demangledSymbol[0] == '\0')
	{
		size_t length = sizeof(demangledSymbol);
		char* val = 0;
		/*DEMANGLER_RESULT result = */demangleToBuffer(symbolName, demangledSymbol, &length, DMGL_MODE_COMPLETE, &val);
		//if(result == DMGL_OK)
	}
}

/**
 * Returns the raw (possibly mangled) symbol name.
 */
template<class Base>
const char* VxSymbolResolver<Base>::mangledSymbol()
{
	resolve();
	return symbolName != 0 ? const_cast<const char*>(symbolName) : 0;
}

/**
 * Returns the demangled symbol name.
 */
template<class Base>
const char* VxSymbolResolver<Base>::symbol()
{
	resolve();
	demangle();
	return demangledSymbol;
}

/**
 * Return the offset in bytes between the function beginning (in machine code)
 * and the stack-frame return address.
 */
template<class Base>
intptr_t VxSymbolResolver<Base>::displacement()
{
	resolve();
	return reinterpret_cast<int>(this->address()) - symbolValue;
}

#endif // defined(VXWORKS)

} /* namespace stack */
} /* namespace utl */
