/* =============================================================================
 * StackWalker.cpp
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

#include "StackWalker.h"

#include <boost/preprocessor/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>

#if defined(linux)
#include <execinfo.h>
#endif

#ifndef MAX_STACKTRACE_DEPTH
#define MAX_STACKTRACE_DEPTH 64
#endif

namespace utl
{
namespace stack
{

#if defined(linux)
size_t GccBackTraceStackWalker::dumpInto(void** begin, void** end)
{
	return backtrace(begin, (end - begin) / sizeof(void*));
}
#endif

/**
 * Preprocessor predicate checking whether the preprocessor loop on
 * stack trace depth has reached the maximum.
 */
#define NOT_MAX_STACK(r, state) \
		BOOST_PP_NOT_EQUAL(state, MAX_STACKTRACE_DEPTH)

/**
 * Preprocessor operator incrementing the stack trace depth counter by one.
 */
#define INC_STACK(r, state) \
		 BOOST_PP_INC(state)

#define ADDR_CLAUSE(func, i) case i: addr = func(i); break;
#define GET_FRAME_ADDRESS(r, state)   ADDR_CLAUSE(__builtin_frame_address, state)
#define GET_RETURN_ADDRESS(r, state)  ADDR_CLAUSE(__builtin_return_address, state)

/**
 * Return the frame address at the specified depth of the call-stack.
 */
void* GccBuiltinStackWalker::getFrameAddress(int i)
{
	void* addr = 0;
	switch(i)
	{
	// calling __builtin_frame_address is only possible with a constant
	// parameter, so we generate a huge ugly switch that handles all cases
	// from 1 to MAX_STACKTRACE_DEPTH.
	BOOST_PP_FOR(0, NOT_MAX_STACK, INC_STACK, GET_FRAME_ADDRESS)
	}
	return addr;
}

/**
 * Return the return address at the specified depth of the call-stack.
 */
void* GccBuiltinStackWalker::getReturnAddress(int i)
{
	void* addr = 0;
	switch(i)
	{
	// calling __builtin_return_address is only possible with a constant
	// parameter, so we generate a huge ugly switch that handles all cases
	// from 1 to MAX_STACKTRACE_DEPTH.
	BOOST_PP_FOR(0, NOT_MAX_STACK, INC_STACK, GET_RETURN_ADDRESS)
	}
	return addr;
}

/**
 * Default constructor.
 */
GccBuiltinStackWalker::iterator::iterator()
	: currentDepth(0), frameAddress(getFrameAddress(currentDepth))
{
}

/**
 * Constructor with a specified frame address.
 */
GccBuiltinStackWalker::iterator::iterator(void* frameAddr)
	: currentDepth(0), frameAddress(frameAddr)
{
}

/**
 * Returns the address of the current position in the call-stack.
 */
void* GccBuiltinStackWalker::iterator::operator*()
{
	return getReturnAddress(currentDepth);
}

/**
 * Go one step deeper into the call-stack.
 */
GccBuiltinStackWalker::iterator& GccBuiltinStackWalker::iterator::operator++()
{
	++currentDepth;
	frameAddress = getFrameAddress(currentDepth);
	return *this;
}

/**
 * Compare iterators. If their frame address is the same, they are considered equal.
 */
bool GccBuiltinStackWalker::iterator::operator==(const iterator& rvalue) const
{
	return frameAddress == rvalue.frameAddress;
}

/**
 * Returns an iterator to the top of the call-stack.
 */
GccBuiltinStackWalker::iterator GccBuiltinStackWalker::begin()
{
	return iterator();
}

/**
 * Returns an iterator representing the bottom of the call-stack.
 */
GccBuiltinStackWalker::iterator GccBuiltinStackWalker::end()
{
	return iterator(0);
}

#if defined(_WIN32)
bool WinStackWalker::isInitialized = false;

void WinStackWalker::libInit()
{
	if(!isInitialized)
	{
		SymSetOptions(SYMOPT_LOAD_LINES);
		SymInitialize(GetCurrentProcess(), NULL, TRUE);
		isInitialized = true;
	}
}

/**
 * Returns an iterator to the top of the call-stack.
 */
WinStackWalker::iterator WinStackWalker::begin()
{
	libInit();
	return iterator(iterator::BEGIN);
}

/**
 * Returns an iterator representing the bottom of the call-stack.
 */
WinStackWalker::iterator WinStackWalker::end()
{
	return iterator(iterator::END);
}

/**
 * Default constructor.
 */
WinStackWalker::iterator::iterator()
{
}

/**
 * Constructor for call-stack beginning or end .
 */
WinStackWalker::iterator::iterator(Type t)
	: context(), stackFrame()
{
	ZeroMemory(&context, sizeof(CONTEXT));
	ZeroMemory(&stackFrame, sizeof(STACKFRAME64));
	if(t == BEGIN)
	{
		next();
	}
}

/**
 * Initializes stack walk.
 */
void WinStackWalker::iterator::beginStackWalk()
{
//#ifdef _M_IX86
//	ZeroMemory(&context, sizeof(CONTEXT));
//	context.ContextFlags = CONTEXT_CONTROL;
//
//	__asm__ __volatile__
//	(
//		"Label:\n"
//		"movl %%ebp, %0\n\t"
//		"movl %%esp, %1\n\t"
//		"movl Label, %%eax\n\t"
//		"movl %%eax, %2\n"
//		: "=m" (context.Ebp), "=m" (context.Esp), "=m" (context.Eip)
//		:
//		: "%eax"
//	);
//#else
	RtlCaptureContext(&context);
//#endif

	stackFrame.AddrPC.Mode       = AddrModeFlat;
	stackFrame.AddrFrame.Mode    = AddrModeFlat;
	stackFrame.AddrStack.Mode    = AddrModeFlat;
#ifdef _M_IX86
	stackFrame.AddrPC.Offset     = context.Eip;
	stackFrame.AddrFrame.Offset  = context.Ebp;
	stackFrame.AddrStack.Offset  = context.Esp;
#elif _M_X64
	stackFrame.AddrPC.Offset     = context.Rip;
	stackFrame.AddrFrame.Offset  = context.Rsp;
	stackFrame.AddrStack.Offset  = context.Rsp;
#elif _M_IA64
	stackFrame.AddrPC.Offset     = context.StIIP;
	stackFrame.AddrFrame.Offset  = context.IntSp;
	stackFrame.AddrBStore.Offset = context.RsBSP;
	stackFrame.AddrBStore.Mode   = AddrModeFlat;
	stackFrame.AddrStack.Offset  = context.IntSp;
#else
#   error "Unsupported platform"
#endif
}

/**
 * Move one step deeper into the call-stack.
 */
void WinStackWalker::iterator::next()
{
	int callerDistance;
	if(stackFrame.AddrPC.Offset == 0)
	{
		callerDistance = 5;
		beginStackWalk();
	}
	else
	{
		callerDistance = 1;
	}
	// the call-stack leading to this function may be differ, whether it is
	// called from the constructor or from the increment operator.
	// When calling this for the first time, we have to go back to the invariant
	// point where iteration takes place, otherwise StackWalk64 will complain
	// the call-stack has changed between calls.
	for(int i = 0; i < callerDistance; ++i)
	{
		WINBOOL result = StackWalk64(
			MACHINE_TYPE,
			GetCurrentProcess(),
			GetCurrentThread(),
			&stackFrame,
			MACHINE_TYPE == IMAGE_FILE_MACHINE_I386 ? 0 : &context,
			0,
			SymFunctionTableAccess64,
			SymGetModuleBase64,
			0);

		if(result == FALSE or stackFrame.AddrPC.Offset == 0)
		{
			stackFrame.AddrPC.Offset = 0;
			break;
		}
	}
}

/**
 * Returns the address of the current position in the call-stack.
 */
void* WinStackWalker::iterator::operator*()
{
	return reinterpret_cast<void*>(stackFrame.AddrPC.Offset);
}

/**
 * Move one step deeper into the call-stack.
 */
WinStackWalker::iterator& WinStackWalker::iterator::operator++()
{
	next();
	return *this;
}

/**
 * Compare iterators. If their PC address is the same, they are considered equal.
 */
bool WinStackWalker::iterator::operator==(const iterator& rvalue) const
{
	return stackFrame.AddrPC.Offset == rvalue.stackFrame.AddrPC.Offset;
}
#endif // defined(_WIN32)

}
}
