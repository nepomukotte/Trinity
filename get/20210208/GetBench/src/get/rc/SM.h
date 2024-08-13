/**
 * @file SM.h
 * @date 9 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors:
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef get_rc_SM_h_INCLUDED
#define get_rc_SM_h_INCLUDED

#include <iostream>

namespace get
{
/// Run Control
namespace rc
{

struct SM
{
	enum State
	{
		OFF     = 0,
		IDLE    = 1,
		DESCRIBED = 2,
		PREPARED  = 3,
		READY   = 4,
		RUNNING = 5,
		PAUSED  = 6
	};

	enum Error
	{
		NO_ERR         = 0,
		WHEN_DESCRIBE  = 1,
		WHEN_PREPARE   = 2,
		WHEN_CONFIGURE = 3,
		WHEN_START     = 4,
		WHEN_STOP      = 5,
		WHEN_PAUSE     = 6,
		WHEN_RESUME    = 7,
		WHEN_BREAKUP   = 8,
		WHEN_RESET	   = 9
	};

	enum Transition
	{
		NO_TRANS        = 0,
		DESCRIBING		= 1,
		PREPARING   	= 2,
		CONFIGURING 	= 3,
		STARTING     	= 4,
		STOPPING     	= 5,
		PAUSING     	= 6,
		RESUMING    	= 7,
		BREAKING_UP   	= 8,
		RESETTING	   	= 9
	};

	struct Status
	{
		State state;
		Error error;
	};

	struct Exception : public Status
	{
		Exception(const SM::Status &, const std::string & msg);
		const std::string errorMessage;
	};
};

std::ostream & operator<<(std::ostream& os, const SM::State & state);
std::ostream & operator<<(std::ostream & os, const SM::Error & error);

} /* namespace rc */
} /* namespace get */
#endif /* get_rc_SM_h_INCLUDED */
