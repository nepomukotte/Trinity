/**
 * @file SM.cpp
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

#include "SM.h"

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
std::ostream & operator<<(std::ostream& os, const SM::State & state)
{
	switch (state)
	{
	case SM::OFF:
		os << "OFF";
		break;
	case SM::IDLE:
		os << "IDLE";
		break;
	case SM::DESCRIBED:
		os << "DESCRIBED";
		break;
	case SM::PREPARED:
		os << "PREPARED";
		break;
	case SM::READY:
		os << "READY";
		break;
	case SM::RUNNING:
		os << "RUNNING";
		break;
	case SM::PAUSED:
		os << "PAUSED";
		break;
	default:
		os << "State(" << int(state) << ')';
		break;
	}
	return os;
}
//_________________________________________________________________________________________________
std::ostream & operator<<(std::ostream & os, const SM::Error & error)
{
	switch (error)
	{
	case SM::NO_ERR:
		os << "NO_ERR";
		break;
	case SM::WHEN_DESCRIBE:
		os << "WHEN_DESCRIBE";
		break;
	case SM::WHEN_PREPARE:
		os << "WHEN_PREPARE";
		break;
	case SM::WHEN_CONFIGURE:
		os << "WHEN_CONFIGURE";
		break;
	case SM::WHEN_START:
		os << "WHEN_START";
		break;
	case SM::WHEN_STOP:
		os << "WHEN_STOP";
		break;
	case SM::WHEN_PAUSE:
		os << "WHEN_PAUSE";
		break;
	case SM::WHEN_RESUME:
		os << "WHEN_RESUME";
		break;
	case SM::WHEN_BREAKUP:
		os << "WHEN_BREAKUP";
		break;
	case SM::WHEN_RESET:
		os << "WHEN_RESET";
		break;
	default:
		os << "Error(" << int(error) << ')';
		break;
	}
	return os;
}
//_________________________________________________________________________________________________
SM::Exception::Exception(const SM::Status & s, const std::string & msg)
	: SM::Status(s), errorMessage(msg)
{
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
