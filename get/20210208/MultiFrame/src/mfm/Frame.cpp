/*
 * $Id: Frame.cpp 1793 2015-09-18 07:47:23Z psizun $
 * @file Frame.cpp
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

#include "Frame.h"
#include "mfm/BlobHeader.h"
#include "mfm/StandardHeader.h"
#include "mfm/Exception.h"
#include "mfm/FrameDictionary.h"
#include <utl/BinIO.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>
namespace io =  boost::iostreams;
using std::ios_base;

namespace mfm
{
//______________________________________________________________________
/**
 * Default constructor.
 * Constructs a minimal frame of the desired kind.
 * @param _kind Kind of frame to create.
 * @param _byteOrder Byte order of the frame to create.
 */
Frame::Frame(FrameKind _kind, utl::Endianness _byteOrder)
	: serializer_(_byteOrder)
{
	// Build header
	if (BLOB == _kind)
	{
		headerPtr_.reset(new BlobHeader(_byteOrder));
	}
	else
	{
		headerPtr_.reset(new StandardHeader(_byteOrder));
		if (BASIC == _kind)
			headerPtr()->set_itemSize_B(1);
	}

	// Resize buffer
	serializer().setCapacity(header().frameSize_B());
	serializer().set_size_B(header().frameSize_B());

	// Set contents of buffer
	updateHeader();
}
//______________________________________________________________________
/**
 * Constructs a frame with the given header.
 * @param _header Frame header to copy and to initialize the frame from.
 * The frame size is set according to the information in the provided header,
 *  all data beyond the said header being initialize to zero.
 */
Frame::Frame(Header const & _header)
	: serializer_(_header.endianness())
{
	// Build header
	Header* h = _header.clone();
	if (_header.isBlobFrame())
	{
		headerPtr_.reset(dynamic_cast<BlobHeader*>(h));
	}
	else
	{
		headerPtr_.reset(dynamic_cast<StandardHeader*>(h));
	}

	// Resize buffer
	serializer().setCapacity(header().frameSize_B());
	serializer().set_size_B(header().frameSize_B());

	// Set contents of buffer
	updateHeader();
}
//______________________________________________________________________
/**
 * Constructs a frame instance from a Serializer.
 * The data are not actually copied.
 * @param _serializer The serializer from which the frame to create will read and write its data.
 */
Frame::Frame(Serializer const & _serializer)
	: serializer_(_serializer)
{
	// Load header from buffer into cache
	loadHeader();
}
//______________________________________________________________________
/**
 * Copy constructor.
 * The copy shares its buffer with the original frame.
 * To create a deep copy, use the Frame::clone() method.
 * @param _frame Original frame to copy.
 */
Frame::Frame(Frame const & _frame)
	: serializer_(_frame.serializer_)
{
	// Load header from buffer into cache
	loadHeader();
}
//______________________________________________________________________
Frame::~Frame()
{
	;
}
//______________________________________________________________________
/**
 * Assignment operator.
 * The copy shares its buffer with the original frame.
 * To create a deep copy, use the Frame::clone() method.
 * @param r Original frame to copy.
 * @return Returns the actual copy.
 */
Frame& Frame::operator=(const Frame & r)
{
	serializer_ = r.serializer();
	// Load header from buffer into cache
	loadHeader();
	return *this;
}
//______________________________________________________________________
/**
 * Creates a clone of this frame by performing a deep copy.
 * @return Returns a deep copy of this instance.
 *
 * Contrary to the copy constructor, the clone does not share its buffer with the original instance.
 */
Frame* Frame::clone() const
{
	// Create a deep copy of the data
	std::auto_ptr< Serializer > s;
	s.reset(serializer_.clone());

	// Create new frame using copied data
	Frame* f = new Frame(*s);

	return f;
}
//______________________________________________________________________
/*
 * (see inherited doc)
 */
Frame const & Frame::frame() const
{
	return *this;
}
//______________________________________________________________________
/*
 * (see inherited doc)
 */
std::string Frame::fieldName(size_t const & offset_B, size_t const & size_B) const
{
	return findFormat().findHeaderField(offset_B, size_B);
}
//______________________________________________________________________
/*
 * (see inherited doc)
 */
void Frame::findBitField(std::string const & fieldName, std::string const & bitFieldName,
		size_t & pos_b, size_t & width_b) const
{
	return findFormat().findHeaderBitField(fieldName, bitFieldName, pos_b, width_b);
}
//______________________________________________________________________
uint32_t Frame::itemCount() const
{
	return header().itemCount();
}
//______________________________________________________________________
uint64_t Frame::itemSize_B(size_t itemIndex) const
{
	if (header().isBasicFrame())
		return header().itemSize_B();
	else if (header().isLayeredFrame())
		return frameSize_B(itemIndex);
	else
		return (header().frameSize_B() - PrimaryHeader::SPEC_SIZE_B);
}
//______________________________________________________________________
/**
 * Returns position of item in Frame buffer, in Bytes.
 * @param itemIndex Index of the item.
 * @return Offset, in Bytes, where the item starts in frame buffer.
 * @throws mfm::ItemNotFound
 */
uint64_t Frame::itemOffset(size_t itemIndex) const
{
	if (itemIndex >= itemCount())
	{
		throw mfm::ItemNotFound(itemIndex, itemCount());
	}
	uint64_t offset = header().headerSize_B();
	if (header().isBasicFrame())
	{
		offset += itemIndex*header().itemSize_B();
	}
	else
	{
		for (size_t i=0; i < itemIndex; ++i)
			offset += itemSize_B(i);
	}
	return offset;
}
//______________________________________________________________________
/**
 * @fn std::auto_ptr< mfm::Header> & Frame::headerPtr()
 * Returns a reference to the auto pointer holding the cached header of this frame.
 * @return A reference to the auto pointer holding the cached header of this frame.
 */
//______________________________________________________________________
/**
 * Reads encoded header and stores its contents into cached header.
 */
void Frame::loadHeader()
{
	std::auto_ptr<Header> headerPtr = Header::decodeHeader(serializer().inputStream());
	headerPtr_ = headerPtr;
}
//______________________________________________________________________
/**
 * Encodes cached header into serialized buffer.
 */
void Frame::updateHeader()
{
	header().encode(serializer().outputStream());
}
//______________________________________________________________________
/**
 * @fn Serializer & Frame::serializer()
 * Non-const accessor to the buffer of this frame.
 * @return A reference to the buffer of this frame.
 */
//______________________________________________________________________
/**
 * @fn Serializer const & Frame::serializer() const
 * Const accessor to the buffer of this frame.
 * @return A const reference to the buffer of this frame.
 */
//______________________________________________________________________
/**
 * Returns an accessor to an item within this frame with given index.
 * @param itemIndex Index of the desired item.
 * @return Returns the item manipulator.
 */
Item Frame::itemAt(size_t const itemIndex)
{
	if (header().isLayeredFrame())
		throw Exception("Operation itemAt() is not supported for frames of the Layered kind!");
	if (itemIndex >= header().itemCount())
		throw Exception("Accessing non-existing item!");
	return Item(this, itemIndex);
}
//______________________________________________________________________
/**
 * Adds a basic item at the end of the frame.
 * Allocates space for a new item at the end of this basic frame and updates header accordingly.
 * FIXME: This method does not take into account the data reserve yet.
 * @return Returns an accessor to the item added to this frame.
 * @throws Throws an exception for non-basic frames.
 */
Item Frame::addItem()
{
	addItems(1);
	return itemAt(header().itemCount()-1);
}
//______________________________________________________________________
/**
 * Adds n basic items at the end of the frame.
 * Allocates space for a new item at the end of this basic frame and updates header accordingly.
 * FIXME: This method does not take into account the data reserve yet.
 * @param n Number of items to add.
 * @throws Throws an exception for non-basic frames.
 */
void Frame::addItems(const size_t & n)
{
	if (not header().isBasicFrame())
		throw Exception("Operation addItem() only supported for frames of the Basic kind!");

	// Update header
	headerPtr()->addItems(n);

	// Allocate memory for new item
	serializer().setCapacity(header().frameSize_B());
	serializer().set_size_B(header().frameSize_B());

	// Re-encode header into buffer
	updateHeader();

	return;
}
//______________________________________________________________________
/**
 * Computes the offset where the embedded frame with given index starts within the frame.
 * @param frameIndex Index of the embedded frame.
 * @return Returns the offset [Bytes].
 */
uint64_t Frame::frameOffset(size_t const & frameIndex) const
{
	if (not header().isLayeredFrame())
		throw Exception("Operation frameAt() is only supported for frames of the Layered kind!");
	if (frameIndex >= frameCount())
		throw OutOfRangeError();

	// First embedded frame starts after end of header
	uint64_t currentOffset = header().headerSize_B();
	size_t currentIndex = 0;
	std::auto_ptr<PrimaryHeader> primaryHeader;

	// Loop over embedded frames
	while (currentIndex < frameIndex)
	{
		// Decode primary header to find out frame size
		// Hack for VxWorks which does not support assignment from return
		std::auto_ptr<PrimaryHeader> tmp = PrimaryHeader::decodePrimaryHeader(serializer().inputStream(currentOffset));
		primaryHeader = tmp;
		currentOffset += primaryHeader->frameSize_B();
		currentIndex++;
	}
	return currentOffset;
}
//______________________________________________________________________
/**
 * @fn uint32_t Frame::frameCount() const
 * Returns the number of frames directly embedded within this frame (without recursivity).
 * @return Number of sub-frames.
 */
//______________________________________________________________________
/**
 * Returns the size of the embedded frame item with given offset.
 * @param frameIndex Index of the variable size item embedded within this frame.
 * @param frameOffset_B Offset [Bytes] where the embedded frame starts within this frame.
 * @return Returns the size [Bytes] of said embedded  frame.
 */
uint64_t Frame::frameSize_B(size_t const & /* frameIndex */, size_t const & frameOffset_B) const
{
	// Decode frame primary header and return size
	std::auto_ptr<PrimaryHeader> primaryHeader =
			PrimaryHeader::decodePrimaryHeader(serializer().inputStream(frameOffset_B));
	return primaryHeader->frameSize_B();
}
//______________________________________________________________________
/**
 * Returns the size of the embedded frame item with given index.
 * @param frameIndex Index of the variable size item embedded within this frame.
 * @return Returns the size [Bytes] of said embedded  frame.
 */
uint64_t Frame::frameSize_B(size_t const & frameIndex) const
{
	// Get frame offset
	uint64_t offset = frameOffset(frameIndex);
	// Decode frame primary header and return size
	return frameSize_B(frameIndex, offset);
}
//______________________________________________________________________
/** Returns a copy of the frame with given index embedded within this layered frame.
 * @param frameIndex Index of the embedded frame to read.
 * @return Returns a hard copy of the given frame.
 * @throws Throws an exception for non-layered frames or if index is out of range.
 */
std::auto_ptr<Frame> Frame::frameAt(size_t const frameIndex)
{
	// Find offset where embedded frame of interest starts
	uint64_t const offset = frameOffset(frameIndex);
	return Frame::read(serializer().inputStream(offset));
}
//______________________________________________________________________
/**
 * Embeds a frame at the end of the frame.
 * Allocates space for a new embedded frame at the end of this layered frame and updates header accordingly.
 * FIXME: This method does not take into account the data reserve yet.
 * @param embeddedFrame The sub-frame to embed into this frame.
 * @throws Throws an exception for non-layered frames.
 */
void Frame::addFrame(Frame const & embeddedFrame)
{
	if (not header().isLayeredFrame())
		throw Exception("Operation addItem(embeddedFrame) only supported for frames of the Layered kind!");

	uint64_t const newItemSize_B = embeddedFrame.header().frameSize_B();

	// Update header
	headerPtr()->addItem(newItemSize_B);

	// Allocate memory for new item
	serializer().setCapacity(header().frameSize_B());
	serializer().set_size_B(header().frameSize_B());

	// Re-encode header into buffer
	updateHeader();

	// Set contents of new embedded frame
	uint64_t const newItemOffset_B = frameOffset(frameCount()-1);
	Serializer embeddedFrameSerializer(serializer(), newItemSize_B, newItemOffset_B);
	embeddedFrameSerializer.read(embeddedFrame.serializer().inputStream(), newItemSize_B);
}
//______________________________________________________________________
/**
 * Returns field of given size starting at given position.
 * This can be an Item field or a Header field, depending on the field container provided.
 * @param container Parent entity of the field.
 * @param size Size of field, in Bytes.
 * @param pos Offset of field, in Bytes, w.r.t beginning of frame.
 * @return Returns the field found.
 */
Field Frame::field(AbstractFieldContainer const* container, const size_t & pos, const size_t & size)
{
	// TODO: add controls
	return Field(container, serializer(), pos, size);
}
//______________________________________________________________________
/**
 * Returns header field of given size starting at given position.
 * @param size Size of field, in Bytes.
 * @param pos Offset of field, in Bytes, w.r.t beginning of frame.
 * @return Returns the header field.
 */
Field Frame::headerField(const size_t &  pos, const size_t & size)
{
	return Field(this, serializer(), pos, size);
}
//______________________________________________________________________
/**
 * Returns header field of given name (in header).
 * @param name Name of field.
 * @return Field found.
 */
Field Frame::headerField(std::string const & name)
{
	size_t fieldPos, fieldSize;
	findHeaderField(name, fieldPos, fieldSize);
	return headerField(fieldPos, fieldSize);
}
//______________________________________________________________________
FrameDictionary & Frame::dictionary()
{
	return FrameDictionary::instance();
}
//______________________________________________________________________
/** Searches the frame dictionary for the format corresponding to this frame.
 *  If the exact revision format is not found,
 *   returns the latest available revision with the same frame type.
 * @return Returns a frame format with the same frame type as this frame
 * @throws mfm::FormatNotFound
 */
FrameFormat const & Frame::findFormat() const
{
	try
	{
		return dictionary().findFormat(header().frameType(), header().revision());
	}
	catch (const FormatRevisionNotFound & e)
	{
		LOG_WARN() << "Could not find description for revision '" << (short) header().revision() << "' of format '" << header().frameType() << '\'';
		return dictionary().findLatestFormat(header().frameType());
	}
}
//______________________________________________________________________
/**
 * Searches for field with given name within  the frame item description.
 * @param name Name of the field to search for.
 * @param[in] itemIndex Index of the item.
 * @param[out] pos Offset of the field, in Bytes,  w.r.t to the beginning of the frame.
 * @param[out] size Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
void Frame::findItemField(const std::string & name, size_t itemIndex, size_t & pos, size_t & size) const
{
	// Exclude non-basic frame metatypes
	if (header().isLayeredFrame())
		throw mfm::Exception("Operation Frame::findItemField is not supported for layered frames!");

	// Find frame format description
	FrameFormat const & format = findFormat();

	// Search field in format
	format.findItemField(name, pos, size);

	// Add offset from header and previous items
	pos += itemOffset(itemIndex);
}
//______________________________________________________________________
/**
 * Searches for field with given name within  the header description.
 * @param name Name of the field to search for.
 * @param[out] pos Offset of the field, in Bytes,  w.r.t to the beginning of the frame.
 * @param[out] size Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
void Frame::findHeaderField(const std::string & name, size_t & pos, size_t & size) const
{
	// Exclude non-basic frame metatypes
	if (header().isBlobFrame())
		throw mfm::Exception("Operation Frame::findHeaderField is not supported for blob frames!");

	// Find frame format description
	FrameFormat const & format = findFormat();

	// Search field in format
	format.findHeaderField(name, pos, size);
}
//______________________________________________________________________
/**
 * Writes this entire frame to an output stream.
 * @param out Output stream to write to.
 */
void Frame::write(std::ostream & out)
{
	serializer().write(out);
}
//______________________________________________________________________
/**
 * Returns a pointer to the underlying buffer.
 * @param out Output stream to write to.
 */
const Byte* Frame::data() const
{
	return serializer().begin();
}
//______________________________________________________________________
/**
 * Reads from an input stream the number of Bytes corresponding to the header reserve,
 *  e.g. to the non-mandatory section of the header, into the buffer of this Frame.
 * @param dataIn Input stream to read from.
 */
void Frame::readHeaderReserve(std::istream & dataIn)
{
	// Make sure that input stream is in exception mode
	dataIn.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);

	size_t const n = header().headerReserveSize_B();
	Serializer reserveSerializer(serializer(), n, header().headerSize_B() - n);
	reserveSerializer.read(dataIn, n);
}
//______________________________________________________________________
/**
 * Reads from an input stream the number of Bytes corresponding to the data section,
 *  into the buffer of this Frame.
 * @param dataIn Input stream to read from.
 */
void Frame::readData(std::istream & dataIn)
{
	// Make sure that input stream is in exception mode
	dataIn.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);

	size_t const n = header().dataSize_B();
	Serializer dataSerializer(serializer(), n, header().headerSize_B());
	dataSerializer.read(dataIn, n);
}
//______________________________________________________________________
/**
 * Reads an entire Frame from an input stream and returns the Frame created.
 * @param dataIn Input stream to read from.
 * @return Returns the frame created and filled.
 *
 * The mandatory header is first decoded to get the kind of frame and its size.
 * The remaining parts are than copied into the frame created.
 */
std::auto_ptr<Frame> Frame::read(std::istream & dataIn)
{
	// Decode header to find out kind of frame (Basic, Layered or Blob)
	std::auto_ptr<Header> headerPtr = Header::decodeHeader(dataIn);

	// Create frame
	std::auto_ptr<Frame> framePtr(new Frame(*headerPtr.get()));

	// Read header reserve and data section
	framePtr->readHeaderReserve(dataIn);
	framePtr->readData(dataIn);

	return framePtr;
}
//______________________________________________________________________
/**
 * Reads an entire Frame from a char buffer and returns the Frame created.
 * @param begin Beginning of the buffer.
 * @param end End of the buffer.
 * @return Returns the frame created and filled.
 *
 * @note Only the first frame in the buffer is returned. The data are copied into the returned frame.
 */
std::auto_ptr<Frame> Frame::read(const char* begin, const char* end)
{
	// Construct stream from char array
	io::array_source source(begin, end);
	io::stream_buffer< io::array_source > streambuf(source);
	std::istream datastream(&streambuf);

	// Read first frame from stream
	std::auto_ptr<Frame> framePtr = mfm::Frame::read(datastream);

	return framePtr;
}
//______________________________________________________________________
/**
 * Skips a given number of successive frames.
 * @param stream Input stream.
 * @param n Number of successive frame(s) to skip.
 * @return Reference to modified stream.
 */
std::istream & Frame::seekFrame(std::istream & stream, const size_t n)
{
	for (size_t i=0; i < n; ++i)
	{
		std::auto_ptr< PrimaryHeader > header = PrimaryHeader::decodePrimaryHeader(stream);
		size_t skipSize_B = header->frameSize_B() - PrimaryHeader::SPEC_SIZE_B;
		stream.seekg(skipSize_B, std::ios_base::cur);
	}
	return stream;
}
//______________________________________________________________________
/**
 * Creates a minimal frame (e.g. without items) corresponding to the given format.
 * @param format The FrameFormat to use.
 * @return Returns a new Frame respecting this format.
 */
std::auto_ptr<Frame> Frame::create(FrameFormat const & format)
{
	// Create header to find out kind of frame (Basic, Layered or Blob)
	std::auto_ptr<Header> headerPtr = Header::createHeader(format);

	// Create frame
	std::auto_ptr<Frame> framePtr(new Frame(*headerPtr.get()));

	return framePtr;
}
//______________________________________________________________________
} /* namespace mfm */
