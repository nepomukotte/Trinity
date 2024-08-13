#include "get/viewgraph.h"
#include "get/plot.h"

#include <QRegExp>
#include <QToolBar>
#include <QToolButton>
#include <QLabel>
#include <QLayout>
#include <QStatusBar>
#include <QPrinter>
#include <QPicture>
#include <QPainter>
#include <QFileDialog>
#include <QImageWriter>
#include <QPrintDialog>
#include <QFileInfo>
#include <QAction>
#include <QKeySequence>

#include <qwt_global.h>
#if _QWT_VERSION < 0x060000
#ifdef QT_SVG_LIB
#include <QSvgGenerator>
#endif
#endif

#include <qwt_counter.h>
#include <qwt_picker_machine.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#if QWT_VERSION >= 0x060000
#include <qwt_plot_renderer.h>
#endif
#include <qwt_plot_curve.h>
#include <qwt_plot_canvas.h>
#include <qwt_text.h>
#include <qwt_math.h>

#include <boost/algorithm/string/replace.hpp>
namespace ba = boost::algorithm;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

class Zoomer: public QwtPlotZoomer
{
public:
	Zoomer(int xAxis, int yAxis, QwtPlotCanvas *canvas) :
			QwtPlotZoomer(xAxis, yAxis, canvas)
	{
		setTrackerMode(QwtPicker::AlwaysOff);
		setRubberBand(QwtPicker::NoRubberBand);

		// RightButton: zoom out by 1
		// Ctrl+RightButton: zoom out to full size

		setMousePattern(QwtEventPattern::MouseSelect2, Qt::RightButton, Qt::ControlModifier);
		setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);
	}
};

#ifndef QT_NO_PRINTER
void ViewGraph::print()
{
	QPrinter printer(QPrinter::HighResolution);

	QString docName = d_plot->title().text();
	if (!docName.isEmpty())
	{
		docName.replace(QRegExp(QString::fromLatin1("\n")), tr(" -- "));
		printer.setDocName(docName);
	}

	printer.setCreator("ACTAR-GET/GetController");
	printer.setOrientation(QPrinter::Landscape);

	QPrintDialog dialog(&printer);
	if (dialog.exec())
	{
#if QWT_VERSION >= 0x060000
		QwtPlotRenderer renderer;

		if (printer.colorMode() == QPrinter::GrayScale)
		{
			renderer.setDiscardFlag(QwtPlotRenderer::DiscardCanvasBackground);
			renderer.setLayoutFlag(QwtPlotRenderer::FrameWithScales);
		}

		renderer.renderTo(d_plot, printer);
#else
		d_plot->print(printer);
#endif
	}
}
#endif // QT_NO_PRINTER

/**
 * Constructs a file basename for export from the title.
 */
QString ViewGraph::buidBaseName() const
{
	std::string basename = d_plot->title().text().toStdString();
	ba::replace_all(basename, " ", "_");
	ba::replace_all(basename, "/", "_");
	ba::replace_all(basename, ".", "_");
	if (basename.empty()) basename = "cobo_frame_display";
	return QString::fromStdString(basename);
}

bool ViewGraph::exportDocument(const QString & fileName)
{
	if (!fileName.isEmpty())
	{
#if QWT_VERSION >= 0x060000

		QwtPlotRenderer renderer;
		renderer.renderDocument(d_plot, fileName, QSizeF(300, 200), 85);

#else // QWT_VERSION <= 5
#ifdef QT_SVG_LIB

		QSvgGenerator generator;
		generator.setFileName(fileName);
		generator.setSize(QSize(800, 600));
		d_plot->print(generator);

#endif // QT_SVG_LIB
#endif // QWT_VERSION
		return true;
	}
	return false;
}

void ViewGraph::exportDocument()
{
	QString baseName = buidBaseName();

#if QWT_VERSION >= 0x060000
	QString suffix(".png");
	QString selFilter;
	QString fileName = baseName + suffix;
#ifndef QT_NO_FILEDIALOG
	const QList<QByteArray> imageFormats = QImageWriter::supportedImageFormats();

	QStringList filters;
	filters += "PDF Documents (*.pdf)";
#ifndef QWT_NO_SVG
	filters += "SVG Documents (*.svg)";
#endif // QWT_NO_SVG
	filters += "Postscript Documents (*.ps)";

	if (imageFormats.size() > 0)
	{
		for (int i = 0; i < imageFormats.size(); i++)
		{
			const QString format(imageFormats[i]);
			const QString filter = QString("%1 Images (*.%2)").arg(format.toUpper()).arg(format);
			if (format == "png") selFilter = filter;
			filters += filter;
		}
	}

	fileName = QFileDialog::getSaveFileName(this, "Export File Name", fileName, filters.join(";;"), &selFilter);
#endif // QT_NO_FILEDIALOG
	exportDocument(fileName);

#else // QWT_VERSION <= 5
#ifdef QT_SVG_LIB
	QString fileName = QString("%1.svg").arg(baseName);
#ifndef QT_NO_FILEDIALOG
	fileName = QFileDialog::getSaveFileName(this, "Export File Name", fileName, tr("SVG images (%1)").arg("*.svg"));
#endif // QT_NO_FILEDIALOG
	exportDocument(fileName);
#endif // QT_SVG_LIB
#endif // QWT_VERSION
}

void ViewGraph::enableZoomMode(bool on)
{
	if (on)
	{
		d_panner->setEnabled(on);

		d_zoomer[0]->setEnabled(on);
		d_zoomer[0]->setZoomBase(true);
		d_zoomer[0]->zoom(0);

		d_zoomer[1]->setEnabled(on);
		d_zoomer[0]->setZoomBase(false);
		d_zoomer[1]->zoom(0);

		d_picker->setEnabled(!on);
	}
	else
	{
		d_plot->setAxisAutoScale(QwtPlot::yLeft);
		d_plot->setAxisAutoScale(QwtPlot::xBottom);

	}
	showInfo();
}

void ViewGraph::zoomOut()
{
	enableZoomMode(false);
}

void ViewGraph::showInfo(QString text)
{
	if (text == QString::null)
	{
		if (d_picker->rubberBand())
			text = "Cursor Pos: Press left mouse button in plot region";
		else
			text = "Zoom: Press mouse button and drag";
	}
	emit statusMessage(text);
}

void ViewGraph::moved(const QPoint &pos)
{
	QString info;
	info.sprintf("Freq=%g, Ampl=%g, Phase=%g", d_plot->invTransform(QwtPlot::xBottom, pos.x()),
			d_plot->invTransform(QwtPlot::yLeft, pos.y()), d_plot->invTransform(QwtPlot::yRight, pos.y()));
	showInfo(info);
}

void ViewGraph::selected()
{
	showInfo();
}

void ViewGraph::setAxisTitle(int axisId, QString title)
{
	d_plot->setAxisTitle(axisId, title);
}

ViewGraph::ViewGraph(QWidget *parent) :
		QWidget(parent)
{
	d_plot = new Plot(tr("Graph"), this);
	setContextMenuPolicy(Qt::NoContextMenu);

	QwtPlotCanvas* canvas = qobject_cast< QwtPlotCanvas* >(d_plot->canvas());

	d_zoomer[0] = new Zoomer(QwtPlot::xBottom, QwtPlot::yLeft, canvas);
	d_zoomer[0]->setRubberBand(QwtPicker::RectRubberBand);
	d_zoomer[0]->setRubberBandPen(QColor(Qt::green));
	d_zoomer[0]->setTrackerMode(QwtPicker::ActiveOnly);
	d_zoomer[0]->setTrackerPen(QColor(Qt::white));

	d_zoomer[1] = new Zoomer(QwtPlot::xTop, QwtPlot::yRight, canvas);

	d_panner = new QwtPlotPanner(d_plot->canvas());
	d_panner->setMouseButton(Qt::MidButton);

#if QWT_VERSION >= 0x060000
	d_picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft, QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn, d_plot->canvas());
	d_picker->setStateMachine(new QwtPickerDragPointMachine());
#else
	d_picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft, QwtPicker::PointSelection | QwtPicker::DragSelection, QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn, d_plot->canvas());
#endif
	d_picker->setRubberBandPen(QColor(Qt::green));
	d_picker->setRubberBand(QwtPicker::CrossRubberBand);
	d_picker->setTrackerPen(QColor(Qt::white));

	QToolBar *toolBar = new QToolBar(this);

	QAction* btnZoom = new QAction(tr("Zoom Out"), toolBar);
	btnZoom->setShortcuts(QKeySequence::ZoomOut);
#if QT_VERSION >= 0x040600
	btnZoom->setIcon(QIcon::fromTheme("zoom-out"));
#endif
	toolBar->addAction(btnZoom);
	connect(btnZoom, SIGNAL(triggered()), SLOT(zoomOut()));

	QAction* actionToggle = new QAction(tr("Toggle All Curves"), this);
	actionToggle->setShortcuts(QKeySequence::SelectAll);
#if QT_VERSION >= 0x040600
	actionToggle->setIcon(QIcon::fromTheme("edit-select-all"));
#endif
	connect(actionToggle, SIGNAL(triggered()), this, SLOT(toggleAllCurves()));
	toolBar->addAction(actionToggle);

	toolBar->addSeparator();

#ifndef QT_NO_PRINTER
	QAction *btnPrint = new QAction(tr("Print"), toolBar);
	btnPrint->setShortcuts(QKeySequence::Print);
#if QT_VERSION >= 0x040600
	btnPrint->setIcon(QIcon::fromTheme("document-print"));
#endif
	toolBar->addAction(btnPrint);
	connect(btnPrint, SIGNAL(triggered()), SLOT(print()));
#endif

	QAction* btnExport = new QAction(tr("Save As"), toolBar);
	btnExport->setShortcuts(QKeySequence::SaveAs);
#if QWT_VERSION < 0x060000
#ifndef QT_SVG_LIB
	btnExport->setDisabled(true);
#endif
#endif
#if QT_VERSION >= 0x040600
	btnExport->setIcon(QIcon::fromTheme("document-save"));
#endif
	toolBar->addAction(btnExport);
	connect(btnExport, SIGNAL(triggered()), SLOT(exportDocument()));

	enableZoomMode(false);
	showInfo();

	d_plot->setAutoReplot(true);
	d_plot->replot();

	connect(d_picker, SIGNAL(moved(const QPoint &)), SLOT(moved(const QPoint &)));
#if QWT_VERSION >= 0x060000
	connect(d_picker, SIGNAL(selected(const QPolygon &)), SLOT(selected()));
#else
	connect(d_picker, SIGNAL(selected(const QwtPolygon &)), SLOT(selected()));
#endif

	QBoxLayout* layout = new QVBoxLayout;
	layout->addWidget(toolBar);
	layout->addWidget(d_plot);
	setLayout(layout);
	setWindowTitle(tr("Graphs"));
}

ViewGraph::~ViewGraph()
{
	clear();
}

void ViewGraph::clear()
{
	d_plot->clear();
}

void ViewGraph::toggleAllCurves()
{
	d_plot->toggleAllCurves();
}

void ViewGraph::updateAxes()
{
#if QWT_VERSION >= 0x060000
	d_plot->updateAxes();
#endif
}

void ViewGraph::setTitle(QString title)
{
	d_plot->setTitle(title);
}

QString ViewGraph::getTitle() const
{
	return d_plot->title().text();
}

int ViewGraph::nbCurve()
{
	int val = d_plot->size();
	return val;
}

void ViewGraph::addCurve(const QString & title, const QString & symbol, const bool yAxisPos, const int offset,
		QVector<double> & xData, QVector<double> & yData)
{
	d_plot->addCurve(title, symbol, yAxisPos, offset, xData, yData);
	updateAxes();
}

void ViewGraph::addCurve(const QString & title, const QString & symbol, const bool yAxisPos, const int size, const double* xData, const double* yData)
{
	d_plot->addCurve(title, symbol, yAxisPos, size, xData, yData);
	updateAxes();
}

void ViewGraph::validatePlot()
{
	d_plot->setAutoReplot(true);
	d_plot->replot();
}
