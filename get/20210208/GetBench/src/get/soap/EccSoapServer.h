/**
 * @file EccSoapServer.h
 * @date Jan 13, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: EccSoapServer.h 1259 2014-04-25 12:17:01Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetBench software project.
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

#ifndef get_soap_EccServer_h_INCLUDED
#define get_soap_EccServer_h_INCLUDED

#include "get/rc/BackEnd.h"
#include "get/rc/SM.h"
#include "ecceccObject.h"
#include ICEUTIL_HEADER(RecMutex.h)
#include <string>

namespace get
{
namespace soap
{
//__________________________________________________________________________________________________
class EccServer
{
public:
	EccServer();
	virtual ~EccServer() {}
	void setConfigManager(const std::string & url=".", const std::string & user="", const std::string & pwd="");
	virtual void createAlarmLogger(const ::utl::net::SocketAddress & address = ::utl::net::SocketAddress(Default::alarmLoggerPortNum));
	void setLoggingBackend();
	/// @name ECC/RCC interface
	/// @{
	int getState(ecc__ResponseState & rs);
	int describe(std::string xmlConfigId, std::string /* options */, ecc__ResponseText & response);
	int prepare(std::string xmlConfigId, std::string /* options */, ecc__ResponseText & response);
	int configure(std::string xmlConfigId, std::string /* options */, ecc__ResponseText & response);
	int specialCmd(std::string command, std::string /* options */, ecc__ResponseText & rt);
	int getConfigIDs(ecc__ResponseText & respText);
	int undo(ecc__ResponseText & response);
	int start(ecc__ResponseText & response);
	int stop(ecc__ResponseText & response);
	int pause(ecc__ResponseText & response);
	int resume(ecc__ResponseText & response);
	int breakup(ecc__ResponseText & response);
	/// @}
protected:
	void updateCachedState();
	ecc__ResponseState cachedState() const;
	bool isBusy() const;
	bool setBusy(const rc::SM::Transition & t, int & code, std::string & msg);
private:
	rc::BackEnd backEnd_;
	ecc__ResponseState cachedState_; ///< Last value of state machine
	IceUtil::RecMutex statusMutex_;
};
//__________________________________________________________________________________________________
} /* namespace soap */
} /* namespace get */
#endif /* get_soap_EccServer_h_INCLUDED */
