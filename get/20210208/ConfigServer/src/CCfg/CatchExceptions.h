/*
 * Except.h
 *
 *  Created on: 6 juil. 2009
 *      Author: rvloic
 */

#ifndef EXCEPT_H_
#define EXCEPT_H_

#include "CCfg/Exception.h"

#if 0 // TODO: replace by a "ICE_BUILD" macro
#include "ConfigServer.h"
#define RETHROW_AS_ICE_ERROR \
	catch(CCfg::Exception& e) \
	{ \
		ExceptionForIce	exIce;		\
		exIce.reason = e.what();	\
		throw (exIce); 				\
	}
#else
// TODO: redefine the macro for SWIG
#define RETHROW_AS_ICE_ERROR		\
	catch(...) 		\
	{								\
		throw;						\
	}

#endif
#define HANDLE_EXCEP(instruction) \
	try \
	{ \
		instruction; \
	} \
	RETHROW_AS_ICE_ERROR
#endif /* EXCEPT_H_ */
