/* =====================================================================================================================
 * CatchAllExceptions.h
 * ---------------------------------------------------------------------------------------------------------------------
 * Macros for catching all exceptions.
 * Created on: Jan 19, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GetController package, whose purpose is to ...
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


#ifndef CatchAllExceptions_h_INCLUDED
#define CatchAllExceptions_h_INCLUDED

#include "mdaq/utl/CmdException.h"
#include <Ice/Ice.h>
#include <exception>

///**********************************************************************
//Definit les exceptions courantes vers ECCServer                    ****
//***********************************************************************
#define CATCH_ALL_EXCEPTIONS()       	                                                \
	catch (const Ice::ConnectionRefusedException &)                                     \
    {                                                                                   \
		LOG_ERROR() << "Connection refused: try restarting the ECC server or check the endpoints"; \
	}                                                                                   \
	catch (const ::mdaq::utl::CmdException & e)                                         \
    {                                                                                   \
		LOG_ERROR() << e.reason;                                                        \
	}                                                                                   \
	catch (const Ice::Exception & e)                                                    \
    {                                                                                   \
		LOG_ERROR() << e.what();                                                        \
	}                                                                                   \
	catch (const std::exception& e)                                                     \
    {                                                                                   \
		LOG_ERROR() << e.what();                                                        \
	}                                                                                   \
	catch (const char* msg)                                                             \
    {                                                                                   \
		LOG_ERROR() << msg;                                                             \
	}                                                                                   \
	catch (const std::string& msg)                                                      \
    {                                                                                   \
		LOG_ERROR() << msg;                                                             \
	}                                                                                   \
	catch (...)                                                                         \
    {                                                                                   \
		LOG_ERROR() << "Unknown exception";                                             \
	}

#endif /* CatchAllExceptions_h_INCLUDED */
