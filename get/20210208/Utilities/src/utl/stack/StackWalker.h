/* =============================================================================
 * StackWalker.h
 * -----------------------------------------------------------------------------
 * class utl::stack::GccBackTraceStackWalker
 * class utl::stack::GccBuiltinStackWalker
 * class utl::stack::WinStackWalker
 * class utl::stack::VxStackWalker
 * Created on: 28 août 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef utl_stack_StackWalker_h_INCLUDED
#define utl_stack_StackWalker_h_INCLUDED

#if defined(_WIN32)
#include <windows.h>
#include <dbghelp.h>
#endif

#include <iterator>

namespace utl
{
namespace stack
{

/**
 * This is a base class for stack iterators.
 * It provides a default implementation for operators that can be derived from
 * another one, and inherits the standard helper class.
 */
template<class I>
class StackWalkerIterator : public std::iterator<std::input_iterator_tag, void*>
{
public:
	I operator++(int) { I oldIterator(*this); ++(*this); return oldIterator; }
	bool operator!=(const I& rvalue) const { return !(*static_cast<const I*>(this) == rvalue); }
};

#if defined(linux)
class GccBackTraceStackWalker
{
public:
	size_t dumpInto(void** begin, void** end);
};
#endif

/**
 * Iterate through the call-stack back to the program entry point
 * using gcc built-in functions.
 * It should work on most platforms.
 */
class GccBuiltinStackWalker
{
public:
	class iterator : public StackWalkerIterator<iterator>
	{
	public:
		void* operator*();
		iterator& operator++();
		bool operator==(const iterator& rvalue) const;

	private:
		friend class GccBuiltinStackWalker;

		iterator();
		iterator(void* frameAddr);

		int currentDepth;
		void* frameAddress;
	};

	iterator begin();
	iterator end();

private:
	static void* getReturnAddress(int i);
	static void* getFrameAddress(int i);
};

#if defined(_WIN32)
/**
 * Iterate through the call-stack back to the program entry point
 * using the Windows API.
 */
class WinStackWalker
{
public:
	class iterator : public StackWalkerIterator<iterator>
	{
	public:
		void* operator*();
		iterator& operator++();
		bool operator==(const iterator& rvalue) const;

	private:
		friend class WinStackWalker;

		enum Type { BEGIN, END };
		iterator();
		iterator(Type t);

		void next();

		void beginStackWalk();

		CONTEXT context;
		STACKFRAME64 stackFrame;
	};

	iterator begin();
	iterator end();

private:
	enum MachineType {
#ifdef _M_IX86
		MACHINE_TYPE = IMAGE_FILE_MACHINE_I386
#elif _M_X64
		MACHINE_TYPE = IMAGE_FILE_MACHINE_AMD64
#elif _M_IA64
		MACHINE_TYPE = IMAGE_FILE_MACHINE_IA64
#endif
	};

	static void libInit();
	static bool isInitialized;
};
#endif

#if defined(VXWORKS)
/**
 * Iterate through the call-stack back to the program entry point
 * using the VxWorks API.
 */
class VxStackWalker
{
public:
	class iterator : public StackWalkerIterator<iterator>
	{
	public:
		void* operator*();
		iterator& operator++();
		bool operator==(const iterator& rvalue) const;

	private:
		friend class VxStackWalker;

		enum Type { BEGIN, END };
		iterator();
		iterator(Type t);

		void next();
	};

	iterator begin();
	iterator end();
};
#endif

}
}

#endif /* utl_stack_StackWalker_h_INCLUDED */
