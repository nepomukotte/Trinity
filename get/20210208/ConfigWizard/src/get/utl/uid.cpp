/**
 * @file uid.cpp
 * @date Jan 28, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: uid.cpp 1105 2014-01-29 08:47:35Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the ConfigWizard software project.
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

#include "get/utl/uid.h"
#include "utl/Logging.h"
#include <boost/lexical_cast.hpp>

namespace get {
namespace utl {
//_________________________________________________________________________________________________
ChannelUID::ChannelUID(const std::string & cobo, size_t asad, size_t aget, size_t channel)
	: coboId(cobo), asadId(asad), agetId(aget), channelId(channel)
{
}
//_________________________________________________________________________________________________
ChannelUID::ChannelUID(const std::string & cobo, const std::string & asad, const std::string & aget, const std::string & channel)
	 : coboId(cobo), asadId(0), agetId(0), channelId(0)
{
	try
	{
		asadId = boost::lexical_cast< size_t >(asad);
		agetId = boost::lexical_cast< size_t >(aget);
		channelId = boost::lexical_cast< size_t >(channel);
	}
	catch (boost::bad_lexical_cast & e)
	{
		LOG_ERROR() << e.what();
	}
}
//_________________________________________________________________________________________________
bool operator< (const ChannelUID & l,const ChannelUID & r)
{
	return (l.coboId < r.coboId or l.asadId < r.asadId or l.agetId < r.agetId or l.channelId < r.channelId);
}
//_________________________________________________________________________________________________
bool operator== (const ChannelUID & l,const ChannelUID & r)
{
	return (l.coboId == r.coboId and l.asadId == r.asadId and l.agetId == r.agetId and l.channelId == r.channelId);
}
//_________________________________________________________________________________________________
std::ostream & operator<< (std::ostream & stream, const ChannelUID & o)
{
	stream << "CoBo[" << o.coboId << "].AsAd[" << o.asadId << "].AGET[" << o.agetId << "].chan[" << o.channelId << ']';
	return stream;
}
//_________________________________________________________________________________________________
} // namespace utl
} // namespace get


