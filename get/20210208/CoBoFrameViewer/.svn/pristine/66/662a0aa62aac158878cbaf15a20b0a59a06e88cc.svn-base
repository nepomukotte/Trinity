/**
 * @file CoBoSurfDisplay.cpp
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

#include "get/CoBoSurfDisplay.h"
#include "get/CoBoEvent.h"
#include "utl/Logging.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QAction>
#include <QKeySequence>
#include <QMessageBox>
#include <qwt3d_function.h>
#include <qwt3d_io.h>

namespace get
{

class CoBoSurf : public Qwt3D::Function
{
public:
	CoBoSurf(Qwt3D::SurfacePlot & pw, const get::CoBoEvent & event, const size_t agetIdx, const get::ChannelMask & mask)
	: Qwt3D::Function(pw), event(event), mask(mask), agetIdx(agetIdx)
	{
	}

	double operator()(double time, double fchannel)
	{
		if (time < 0 or time >= 512 or fchannel < 0 or fchannel >=68) return 0.;
		size_t chanIdx = floor(fchannel);
		try
		{
			if (mask.count(chanIdx) <= 0) return 0.;
			size_t timeIdx = floor(time) - event.readOffset();
			try
			{
				return event.channel(chanIdx, agetIdx).sampleAt(timeIdx);
			}
			catch (std::exception & e)
			{
			}
		}
		catch (const get::CoBoEvent::ChannelNotFound & e)
		{
		}
		return 0.;
	}
private:
	const get::CoBoEvent & event;
	const get::ChannelMask & mask;
	size_t agetIdx;
};

//_________________________________________________________________________________________________
CoBoSurfDisplay::CoBoSurfDisplay(QWidget* parentWidget) : Qwt3D::SurfacePlot(parentWidget)
{
	setTitleFont("Courier",12, QFont::Bold);
	coordinates()->setLabelFont("Courier",11, QFont::Bold);
	coordinates()->axes[Qwt3D::X1].setLabelString(tr("Time Bucket Index"));
	coordinates()->axes[Qwt3D::X2].setLabelString(tr("Time Bucket Index"));
	coordinates()->axes[Qwt3D::X3].setLabelString(tr("Time Bucket Index"));
	coordinates()->axes[Qwt3D::X4].setLabelString(tr("Time Bucket Index"));
	coordinates()->axes[Qwt3D::Y1].setLabelString(tr("Channel"));
	coordinates()->axes[Qwt3D::Y2].setLabelString(tr("Channel"));
	coordinates()->axes[Qwt3D::Y3].setLabelString(tr("Channel"));
	coordinates()->axes[Qwt3D::Y4].setLabelString(tr("Channel"));
	coordinates()->axes[Qwt3D::Z1].setLabelString(tr("Digitized Signal Value"));
	coordinates()->axes[Qwt3D::Z2].setLabelString(tr("Digitized Signal Value"));
	coordinates()->axes[Qwt3D::Z3].setLabelString(tr("Digitized Signal Value"));
	coordinates()->axes[Qwt3D::Z4].setLabelString(tr("Digitized Signal Value"));
    for (unsigned i=0; i!=coordinates()->axes.size(); ++i)
    {
    	coordinates()->axes[i].setMajors(4);
    	coordinates()->axes[i].setMinors(5);
    	//coordinates()->axes[i].setLabelColor(Qwt3D::RGBA(0,0,0.4));
    	coordinates()->axes[i].setLabel(true);
    	coordinates()->axes[i].setAutoScale(true);
    	coordinates()->axes[i].setLabelPosition(Qwt3D::Triple(0, 0, 0), Qwt3D::Center);
    }

	setCoordinateStyle(Qwt3D::BOX);
	setFloorStyle(Qwt3D::NOFLOOR);
	showColorLegend(true);
	legend()->setAutoScale(true);
	setRotation(50, 0, -30);
	setViewportShift(0.05,0);
	setOrtho(false);
	enableLighting(true);

	QAction* actionSave = new QAction(tr("Save As "), this);
	actionSave->setShortcuts(QKeySequence::Save);
	connect(actionSave, SIGNAL(triggered()), this, SLOT(saveAs()));
	addAction(actionSave);
	setContextMenuPolicy(Qt::ActionsContextMenu);
}
//_________________________________________________________________________________________________
void CoBoSurfDisplay::saveAs()
{
	QStringList formats = Qwt3D::IO::outputFormatList();
	QString filter = formats.join(" *.").prepend("*.").replace("_", ".").toLower();
	QString defaultName = QFileInfo(windowTitle()).fileName().replace(" ", "_").append(".png");
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"), defaultName,
		tr("Images (%1);; All files (%2)").arg(filter).arg("*"));
	if (not saveAs(fileName))
	{
		QMessageBox::warning(this, tr("Error"), tr("Error exporting image as file '%1'.").arg(fileName));
	}
}
//_________________________________________________________________________________________________
void CoBoSurfDisplay::setData(const get::CoBoEvent & event, const size_t agetIdx, const get::ChannelMask & mask)
{
	if (event.containsAGetChip(agetIdx))
	{
		setPlotStyle(Qwt3D::FILLED); //FILLEDMESH);
		setMeshLineWidth(0.5);

		CoBoSurf func(*this, event, agetIdx, mask);
		func.setMesh(2*512, 3*68);
		func.setDomain(0., 512., 0., 68.);
		func.setMinZ(0.);
		func.create();

		setScale(1., 5., 0.5);

		updateData();
		updateColorLegend();
		updateGL();
	}
	else
	{
		setPlotStyle(Qwt3D::NOPLOT);
	}
}
//_________________________________________________________________________________________________
void CoBoSurfDisplay::updateColorLegend(int majors, int minors)
{
	legend()->setMajors(majors);
	legend()->setMinors(minors);
	double start, stop;
	coordinates()->axes[Qwt3D::Z1].limits(start,stop);
	legend()->setLimits(start, stop);
}
//_________________________________________________________________________________________________
void CoBoSurfDisplay::enableLighting(bool enable)
{
	if (plotStyle() != Qwt3D::NOPLOT)
	{
		SurfacePlot::enableLighting(enable);
		illuminate(0);
		updateGL();
	}
}
//_________________________________________________________________________________________________
void CoBoSurfDisplay::setLabel(const QString & title)
{
	setWindowTitle(title);
	setTitle(title);
}
//_________________________________________________________________________________________________
QString CoBoSurfDisplay::getLabel() const
{
	return windowTitle();
}
//_________________________________________________________________________________________________
bool CoBoSurfDisplay::saveAs(const QString & fileName)
{
	if (not fileName.isEmpty())
	{
		// Find format
		const QStringList formats = Qwt3D::IO::outputFormatList();
		QString format("png");
		for (int i = 0; i < formats.size(); ++i)
		{
			if (fileName.endsWith(formats.at(i).toLower().replace("_", ".")))
			{
				format = formats.at(i);
				break;
			}
		}
		return Qwt3D::IO::save(this, fileName, format);
	}
	return false;
}
//_________________________________________________________________________________________________
} /* namespace get */
