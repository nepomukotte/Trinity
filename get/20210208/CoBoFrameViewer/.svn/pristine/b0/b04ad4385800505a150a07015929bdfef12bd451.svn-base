/**
 * @file CoBoStandardDisplay.cpp
 * @date 16 nov. 2012
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

#include "get/CoBoStandardDisplay.h"
#include "get/Channel.h"
#include "utl/Logging.h"
#include <QSettings>
#include <QString>
#include <qwt_plot.h>
#include <algorithm> // transform
#include <functional> // plus
#include <QtGlobal>

namespace get
{
//_________________________________________________________________________________________________
CoBoStandardDisplay::CoBoStandardDisplay(QWidget* parentWidget) : ViewGraph(parentWidget)
{
	setTitle("");
	setAxisTitle(QwtPlot::xBottom, tr("Time Bucket Index"));
	setAxisTitle(QwtPlot::yLeft, tr("Digitized Signal Value"));
	symbolStyleName = QSettings().value("symbolStyle", "NoSymbol").toString();
	randomColorMode = QSettings().value("randomColorsEnabled", false).toBool();
}
//_________________________________________________________________________________________________
void CoBoStandardDisplay::setData(const get::CoBoEvent & event, const size_t agetIdx, const get::ChannelMask & mask)
{
	try
	{
		// Set contents of Qwt graph
		clear();

		if (event.containsAGetChip(agetIdx))
		{
			typedef std::set<short> ChannelSet;
			ChannelSet::const_iterator chanIter;
			for (chanIter = mask.begin(); chanIter != mask.end(); ++chanIter)
			{
				try
				{
					short chanIdx = *chanIter;
					std::vector< double > xData(event.channel(chanIdx, agetIdx).buckIndexes().begin(),
							event.channel(chanIdx, agetIdx).buckIndexes().end());
					// Add readout offset
					if (event.readOffset() > 0)
					{
						transform(xData.begin(), xData.end(), xData.begin(), bind2nd(std::plus< double > (), event.readOffset()));
					}
					const std::vector< double > yData(event.channel(chanIdx, agetIdx).sampleValues().begin(),
						event.channel(chanIdx, agetIdx).sampleValues().end());
					//const double* xData = &event.channel(chanIdx, agetIdx).buckIndexes().front();
					//const double* yData = &event.channel(chanIdx, agetIdx).sampleValues().front();
					addCurve(event.channel(chanIdx, agetIdx).sampleCount(), &xData.front(), &yData.front(),
							QObject::tr("Channel %1").arg(chanIdx), symbolStyleName, channelColor(chanIdx));
				}
				catch (const get::CoBoEvent::ChannelNotFound & e)
				{
					//LOG_DEBUG() << e.what();
				}
			}
		}
		updateAxes();
		validatePlot();
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
}
//__________________________________________________________________________________________________
QColor CoBoStandardDisplay::channelColor(const size_t channelId) const
{
	if (randomColorMode or channelId >= MAX_COLORS)
		return randomColor();
	return QColor(predefinedColors[channelId]);
}
//__________________________________________________________________________________________________
/**
 * Creates a random color.
 * @return Random color created.
 */
QColor CoBoStandardDisplay::randomColor()
{
	return QColor(qrand()%255, qrand()%255, qrand()%255);
}
//________________________________________________________________________________________________
const QString CoBoStandardDisplay::predefinedColors[] = { "#bf8525", "#128657", "#fcddda",
        "#f99392", "#049ce9", "#87d065", "#d2b91a", "#198f43", "#da5ece", "#8dd0ea", "#2a1e48",
        "#b780b9", "#df4197", "#837f85", "#489544", "#d5e97f", "#15b9e4", "#a5ba99", "#7e0835",
        "#f4496c", "#5765ef", "#628ee0", "#ceddb2", "#49e45f", "#73d84c", "#9f6c8d", "#6f8ed0",
        "#c89cac", "#675ea4", "#ddb638", "#e28b7f", "#b7b9b5", "#7060f1", "#e86be3", "#74abdf",
        "#66d246", "#a9ae36", "#36dfd7", "#3f3725", "#ae912a", "#92f804", "#e7916d", "#29ceba",
        "#53854f", "#44d47d", "#87ad06", "#f5c4ed", "#22af68", "#5c13e0", "#838759", "#bd0470",
        "#b5fa4d", "#893cf5", "#98a663", "#a3e2a4", "#0089fd", "#8ba37c", "#ab0f9d", "#fc4e79",
        "#b8d0de", "#6f8e10", "#3fb94d", "#283f52", "#d2dded", "#2206cc", "#a13822", "#9e60b5"};
//_________________________________________________________________________________________________
QString CoBoStandardDisplay::getLabel() const
{
	return getTitle();
};
//_________________________________________________________________________________________________
void CoBoStandardDisplay::setLabel(const QString & title)
{
	setTitle(title);
	validatePlot();
};
//_________________________________________________________________________________________________
bool CoBoStandardDisplay::saveAs(const QString & fileName)
{
	return exportDocument(fileName);
}
//_________________________________________________________________________________________________
} /* namespace get */
