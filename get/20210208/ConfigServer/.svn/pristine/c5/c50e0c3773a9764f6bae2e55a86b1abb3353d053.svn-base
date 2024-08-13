/**
 * Declarations For C++
 */
%module CCfg
%{	
	#include "src/CCfg/CConfigFacade.h"
	#include "src/CCfg/DocumentFacade.h"
	#include "src/CCfg/CatchExceptions.h"	
%}

/**
 * Preload of the C++ library
 */
%pragma(java) jniclasscode=%{
	static {
		try { System.loadLibrary("ccfg"); }
		catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load :\n   "
				 				+ e + "\n   " +
				 				"at LoadLibrary in SwigFacade.i : Wrong library name");
      		System.exit(1);
		}
	}
%}

/**
 * Manage Exception
 */
%include exception.i
%include "CatchExceptions.h"
%exception
{
	try
	{ $action }
	catch(CCfg::Exception& e)  			//C++ Exception bring to Java Exception	
	{ SWIG_exception(SWIG_SyntaxError, const_cast<char*>(e.what())); } 
}

/**
 * Manage shared_ptr 
 */
%include boost_shared_ptr.i
SWIG_SHARED_PTR(CConfigFacade, CCfg::CConfigFacade)
SWIG_SHARED_PTR(DocumentFacade, CCfg::DocumentFacade)

/**
 * Manage String (wrapper)
 */
%include std_string.i

/**
 * Manage Vector (wrapper)
 */
%include std_vector.i

namespace std { %template(StringVector) vector<string>; };

/*
 * Declarations For Swig
 */
%include "CConfigFacade.h"
%include "DocumentFacade.h"
