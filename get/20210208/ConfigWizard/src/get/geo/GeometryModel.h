/**
 * @file GeometryModel.h
 * @date Jan 24, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: GeometryModel.h 1105 2014-01-29 08:47:35Z psizun $
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

#ifndef get_geo_GeometryModel_h_INCLUDED
#define get_geo_GeometryModel_h_INCLUDED

#include "get/utl/uid.h"
#include <QAbstractTableModel>
#include <QtGlobal>
/*
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
namespace bmi = boost::multi_index;
using bmi::indexed_by;
using bmi::member;
using bmi::multi_index_container;
using bmi::ordered_unique;
using bmi::random_access;
using bmi::sequenced;
using bmi::tag;
*/
#include <vector>

class TH2F;

namespace get
{
// Forward declaration
namespace cfg
{
	class ChannelConfigModel;
}
namespace geo
{
//_________________________________________________________________________________________________
struct PadChannel
{
	utl::ChannelUID channel;
	utl::PadUID pad;
	inline float x() const { return pad.first; }
	inline float y() const { return pad.second; }
};
//_________________________________________________________________________________________________
/*struct Pad {};
struct Channel {};
typedef multi_index_container
<
	PadChannel,
	indexed_by
	<
		random_access<>,
		ordered_unique< tag<Channel>, BOOST_MULTI_INDEX_MEMBER(PadChannel, ChannelUID, channel) >,
		ordered_unique< tag<Pad>, BOOST_MULTI_INDEX_MEMBER(PadChannel, PadUID, pad) >
  	>
> Mapping;
typedef Mapping::index<Pad>::type MappingPadIdx;
typedef Mapping::index<Channel>::type MappingChannelIdx;
*/
typedef std::vector< PadChannel > Mapping;
//_________________________________________________________________________________________________
/**
 * A class to store the mapping between electronics channels and detector plane pads.
 */
class GeometryModel : public QAbstractTableModel
{
public:
	GeometryModel(QObject* parentWidget=0);
	int rowCount(const QModelIndex & parent = QModelIndex()) const;
	int columnCount(const QModelIndex & parent = QModelIndex()) const;
	Qt::ItemFlags flags(const QModelIndex & index = QModelIndex()) const;
	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
	const PadChannel & mapping(const int row) const;
	QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
	bool setData(const QModelIndex & index, const QVariant & value, int role);
	void addChannel(const PadChannel & padChannel = PadChannel());
	void removeLastChannel();
	void clear();
	float xMin() const;
	float xMax() const;
	float yMin() const;
	float yMax() const;
	TH2F* createHistogram();
	TH2F* createThresholdHistogram(cfg::ChannelConfigModel* );
	TH2F* createLsbThresholdHistogram(cfg::ChannelConfigModel* );
private:
	enum Columns
	{
		CoBoColumn, ///< the index of the column displaying the CoBo identifier
		AsAdColumn,	///< the index of the column displaying the AsAd identifier
		AGETColumn, ///< the index of the column displaying the AGET chip identifier
		ChannelColumn,	///< the index of the column displaying the AGET channel number
		XColumn, ///< the index of the column displaying the pad row
		YColumn, ///< the index of the column displaying the pad column
		COLUMN_COUNT	///< the number of column in the model.
	};
	static const std::string columnNames[COLUMN_COUNT];
	Mapping mapping_;
};
//_________________________________________________________________________________________________
} /* namespace geo */
} /* namespace get */
#endif /* get_geo_GeometryModel_h_INCLUDED */
