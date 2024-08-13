/* =============================================================================
 * $Id: Header.cpp 907 2012-11-05 09:16:44Z psizun $
 * @file Header.cpp
 * @author sizun
 * @date 15 déc. 2011
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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

#include "mfm/Header.h"
#include "mfm/BlobHeader.h"
#include "mfm/StandardHeader.h"
#include "mfm/Exception.h"

namespace mfm {
//______________________________________________________________________
/**
 * Default constructor.
 * Constructs a header.
 * @param _byteOrder Byte order.
 */
Header::Header(utl::Endianness _byteOrder)
	: PrimaryHeader(_byteOrder)
{
    ;
}
//______________________________________________________________________
/**
 * Constructs a header based on a primary header.
 * @param primaryHeader The primary header to copy.
 */
Header::Header(PrimaryHeader const & primaryHeader)
    : PrimaryHeader(primaryHeader)
{
    ;
}
//______________________________________________________________________
/**
 * @fn Header::~Header()
 * Destructor.
 */
//______________________________________________________________________
/**
 * @fn Header* Header::clone() const
 * Clones this header.
 * @return Returns the actual clone.
 */
//______________________________________________________________________
/**
 * Reads from a binary input stream and decodes a header.
 * @param dataIn Input stream to read from.
 * @return Returns a new header.
 *
 * First reads and decodes a primary header. If it turns out this primary header does not identify a blob frame,
 *  reads further to decode the remaining (mandatory) header bytes.
 */
std::auto_ptr<Header> Header::decodeHeader(std::istream & dataIn)
{
	std::auto_ptr<PrimaryHeader> primaryHeaderPtr = PrimaryHeader::decodePrimaryHeader(dataIn);

	std::auto_ptr<Header> headerPtr;
	if (primaryHeaderPtr->isBlobFrame())
	{
		headerPtr.reset(new BlobHeader(*primaryHeaderPtr));
	}
	else
	{
		headerPtr.reset(new StandardHeader(*primaryHeaderPtr, dataIn));
	}

	return headerPtr;
}
//______________________________________________________________________
/**
 * Creates a header from scratch appropriate for a frame with given format description.
 * @param format Frame format description.
 * @return Returns a new header.
 */
std::auto_ptr<Header> Header::createHeader(FrameFormat const & format)
{
	std::auto_ptr<PrimaryHeader> primaryHeaderPtr = PrimaryHeader::createPrimaryHeader(format);

	std::auto_ptr<Header> headerPtr;
	if (primaryHeaderPtr->isBlobFrame())
	{
		headerPtr.reset(new BlobHeader(*primaryHeaderPtr));
	}
	else
	{
		headerPtr.reset(new StandardHeader(*primaryHeaderPtr, format));
	}

	return headerPtr;
}
//______________________________________________________________________
/**
 * Updates the header to reflect the addition of an item.
 * Implementation depends on the kind of frame.
 */
void Header::addItem()
{
    throw Exception("Operation addItem() not supported by class Header!");
}
//______________________________________________________________________
/**
 * Updates the header to reflect the addition of n items.
 * Implementation depends on the kind of frame.
 */
void Header::addItems(size_t const & /* n */)
{
    throw Exception("Operation addItems() not supported by class Header!");
}
//______________________________________________________________________
/**
 * Updates the header to reflect the addition of an item of a certain size.
 * Implementation depends on the kind of frame.
 * @param itemSize Size [Bytes] of the item to add.
 */
void Header::addItem(uint64_t const /* itemSize */)
{
    throw Exception("Operation addItem(itemSize) not supported by class Header!");
}
//______________________________________________________________________
} /* namespace mfm */
