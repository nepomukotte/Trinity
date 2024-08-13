/*
 * $Id: Frame.h 1039 2013-02-13 08:34:27Z psizun $
 * @file Frame.h
 * @created 9 févr. 2012
 * @author sizun
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
 * -----------------------------------------------------------------------------
 */

#ifndef mfm_Frame_h_INCLUDED
#define mfm_Frame_h_INCLUDED

#include <mfm/Common.h>
#include <mfm/Serializer.h>
#include <mfm/PrimaryHeader.h>
#include <mfm/Header.h>
#include <mfm/StandardHeader.h>
#include <mfm/BlobHeader.h>
#include <mfm/Item.h>
#include <mfm/Field.h>
#include <mfm/AbstractFieldContainer.h>
#include <memory>
#include <string>

namespace mfm {
//______________________________________________________________________
class FrameDictionary;
class FrameFormat;
//______________________________________________________________________
/**
 * This is the main class of the API. It represents a full frame.
 */
class Frame : public AbstractFieldContainer
{
public:
	Frame(FrameKind _kind =  LAYERED, utl::Endianness _byteOrder = utl::BigEndian);
	Frame(Header const &);
	Frame(Serializer const &);
	Frame(Frame const &);
	virtual ~Frame();
	Frame& operator=(const Frame & r);
    virtual Frame* clone() const;
	Frame const & frame() const;
	Header const & header() const { return *headerPtr_; }
	uint64_t offset_B() const { return 0u; }

	std::string fieldName(size_t const & offset_B, size_t const & size_B) const;
	void findBitField(std::string const & fieldName, std::string const & bitFieldName,
			size_t & pos_b, size_t & width_b) const;

	/** @name Access to items */
	///@{
	uint32_t itemCount() const;
	uint64_t itemSize_B(size_t itemIndex=0) const;
	uint64_t itemOffset(size_t itemIndex) const;
	Item itemAt(size_t const itemIndex);
	Item item() { return itemAt(0); }
	Item addItem();
	void addItems(const size_t & n=1);
	///@}

	/** @name Access to embedded frames
	 * These methods only concern layered frames. */
	///@{
	uint32_t frameCount() const { return itemCount(); }
	uint64_t frameSize_B(size_t const & frameIndex) const;
	uint64_t frameSize_B(size_t const & /* frameIndex */, size_t const & frameOffset_B) const;
	uint64_t frameOffset(size_t const & frameIndex) const;
	std::auto_ptr<Frame> frameAt(size_t const frameIndex);
	void addFrame(Frame const & embeddedFrame);
	///@}

	Field field(AbstractFieldContainer const* , const size_t & pos, const size_t & size);
	Field headerField(const size_t & pos, const size_t & size);
	Field headerField(std::string const & name);

	void readHeaderReserve(std::istream & dataIn);
	void readData(std::istream & dataIn);
	static std::istream & seekFrame(std::istream & in, const size_t n);
	static std::auto_ptr<Frame> read(std::istream & in);
	static std::auto_ptr<Frame> read(const char* begin, const char* end);
	static std::auto_ptr<Frame> create(FrameFormat const & format);
	void write(std::ostream & out);
	const Byte* data() const;
	void findItemField(const std::string & name, size_t itemIndex, size_t & pos, size_t & size) const;
	void findHeaderField(const std::string & name, size_t & pos, size_t & size) const;
	FrameFormat const & findFormat() const;
	Serializer & serializer() { return serializer_; }
private:
	std::auto_ptr< mfm::Header> & headerPtr() { return headerPtr_; }
	void loadHeader();
	void updateHeader();
	Serializer const & serializer() const { return serializer_; }
	/// Returns the dictionary of all known frame formats.
	static FrameDictionary & dictionary();
private:
	Serializer serializer_; ///< Serializer encapsulating all of the frame data.
	std::auto_ptr< mfm::Header> headerPtr_; ///< Cached header
};
//______________________________________________________________________
} /* namespace mfm */
#endif /* mfm_Frame_h_INCLUDED */
