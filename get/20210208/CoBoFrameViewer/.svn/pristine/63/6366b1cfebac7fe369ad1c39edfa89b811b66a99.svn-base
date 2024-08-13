/**
 * @file ChannelMask.cpp
 * @date 26 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the CoBoFrameViewer software project.
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

#include "get/ChannelMask.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
namespace ba = boost::algorithm;
#include <boost/lexical_cast.hpp>
#include <vector>
#include <list>

namespace get
{
//__________________________________________________________________________________________________
ChannelMask::ChannelMask()
{
}
//__________________________________________________________________________________________________
std::string ChannelMask::toString()  const
{
	typedef std::set< short >::const_iterator MaskIter;
	std::list< std::string > tokens;
	for (MaskIter i=begin(); i!= end();)
	{
		short firstChannel = *i, lastChannel=firstChannel;
		for (short nextChannel=firstChannel; nextChannel <= 68; ++nextChannel)
		{
			if (find(nextChannel) != end())
			{
				lastChannel = nextChannel;
			}
			else
			{
				break;
			}
		}
		std::ostringstream token;
		if (firstChannel >= lastChannel)
			token << lastChannel;
		else
			token << firstChannel << '-' << lastChannel;
		tokens.push_back(token.str());
		i = lower_bound(lastChannel+1);
	}
	return ba::join(tokens, ",");
}
//__________________________________________________________________________________________________
ChannelMask ChannelMask::fromString(const std::string & pattern)
{
	ChannelMask mask;
	// Split user input
	std::vector<std::string> tokens;
	ba::split(tokens, pattern, ba::is_any_of(","), ba::token_compress_on);
	for (size_t i = 0; i < tokens.size(); ++i)
	{
		try
		{
			int channel = boost::lexical_cast<int>(tokens[i]);
			mask.insert(channel);
		}
		catch (const boost::bad_lexical_cast & e)
		{
			std::istringstream iss(tokens[i]);
			short first, last;
			char sep;
			if ((iss >> first >> sep >> last) and first <= last)
			{
				for (short i = first; i <= last; ++i)
					mask.insert(i);
			}
		}
	}
	return mask;
}
//__________________________________________________________________________________________________
} /* namespace get */
