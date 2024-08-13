
#include "plot.h"
#include "get/PlotCurve.h"
#include "get/complexnumber.h"

#include <qwt_math.h>
#include <qwt_scale_engine.h>
#include <qwt_symbol.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_global.h>
#if QWT_VERSION < 0x060100
#include <qwt_legend_item.h>
#endif
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <qmath.h>
#include <QVectorIterator>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

#if QT_VERSION < 0x040601
#define qExp(x) std::exp(x)
#define qAtan2(y, x) std::atan2(y, x)
#endif
//__________________________________________________________________________________________________
int Plot::size()
{
	return arrayOfCurves.size();
}
//__________________________________________________________________________________________________
Plot::Plot(QString title, QWidget *parent):
    QwtPlot(parent)
{
	srand ( time(NULL) );
    setAutoReplot(false);
    setTitle(title);
    setCanvasBackground(QColor(Qt::white));

    // Grid
    QwtPlotGrid *grid = new QwtPlotGrid;
    grid->enableXMin(true);
#if QWT_VERSION < 0x060100
    grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->setMinPen(QPen(Qt::gray, 0 , Qt::DotLine));
#else
    grid->setMajorPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->setMinorPen(QPen(Qt::gray, 0 , Qt::DotLine));
#endif
    grid->attach(this);

    // Axes
    enableAxis(QwtPlot::yLeft);
    setAxisTitle(QwtPlot::xBottom, tr("Time"));
    setAxisTitle(QwtPlot::yLeft, tr("Amplitude"));
    setAxisTitle(QwtPlot::yRight, tr("Amplitude 2"));

    setAxisMaxMajor(QwtPlot::xBottom, 6);
    setAxisMaxMinor(QwtPlot::xBottom, 10);
    setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine);

    // Markers
    d_marker1 = new QwtPlotMarker();
    d_marker1->setValue(0.0, 0.0);
    d_marker1->setLineStyle(QwtPlotMarker::VLine);
    d_marker1->setLabelAlignment(Qt::AlignRight | Qt::AlignBottom);
    d_marker1->setLinePen(QPen(Qt::green, 0, Qt::DashDotLine));
    d_marker1->attach(this);

    d_marker2 = new QwtPlotMarker();
    d_marker2->setLineStyle(QwtPlotMarker::HLine);
    d_marker2->setLabelAlignment(Qt::AlignRight | Qt::AlignBottom);
    d_marker2->setLinePen(QPen(QColor(200,150,0), 0, Qt::DashDotLine));
#if QWT_VERSION >= 0x060000
    d_marker2->setSymbol( new QwtSymbol(QwtSymbol::Diamond, QColor(Qt::yellow), QColor(Qt::green), QSize(7,7)));
#else
    d_marker2->setSymbol(QwtSymbol(QwtSymbol::Diamond, QColor(Qt::yellow), QColor(Qt::green), QSize(7,7)));
#endif
    d_marker2->attach(this);

    QwtLegend* legend = new QwtLegend(this);
    insertLegend(legend, QwtPlot::RightLegend, 0.1);
#if QWT_VERSION < 0x060100
    legend->setItemMode(QwtLegend::CheckableItem);
    QObject::connect(this, SIGNAL(legendChecked(QwtPlotItem *, bool)),
        this, SLOT(showCurve(QwtPlotItem *, bool)));
#else
    legend->setDefaultItemMode(QwtLegendData::Checkable);
    QObject::connect(legend, SIGNAL(checked(const QVariant &, bool, int)),
            this, SLOT(showCurve(const QVariant &, bool, int)));
#endif

    setAutoReplot(true);
    replot();
}
//__________________________________________________________________________________________________
Plot::~Plot()
{
	clear();
}
//__________________________________________________________________________________________________
void Plot::clear()
{
	for (int i=0; i<arrayOfCurves.size(); i++)
	{
		arrayOfCurves.at(i)->detach();
		delete arrayOfCurves.at(i);
	}
	arrayOfCurves.clear();
	setTitle("");
}
//__________________________________________________________________________________________________
/*Style {
  NoSymbol = -1,
  Ellipse,
  Rect,
  Diamond,
  Triangle,
  DTriangle,
  UTriangle,
  LTriangle,
  RTriangle,
  Cross,
  XCross,
  HLine,
  VLine,
  Star1,
  Star2,
  Hexagon,
  UserSymbol = 1000
}*/
//__________________________________________________________________________________________________
QwtPlotCurve* Plot::createCurve(const QString & title, const QString & symbolName, const bool ySite)
{
	QwtPlotCurve* curve = new get::PlotCurve(title);
	curve->setTitle(title);
	curve->setRenderHint(QwtPlotItem::RenderAntialiased);

	// Color
	QColor randomColor(randNb(0,255), randNb(0,255), randNb(0,255), 255);

	// Pen
	QPen pen(randomColor);
	curve->setPen(pen);

	// Symbol
	setCurveSymbol(curve, symbolName, pen);

	if (ySite)
		curve->setYAxis(QwtPlot::yLeft);
	else
		curve->setYAxis(QwtPlot::yRight);

	return curve;
}
//__________________________________________________________________________________________________
void Plot::setCurveSymbol(QwtPlotCurve* curve, const QString & symbolName, const QPen & pen)
{
	// Symbol
#if QWT_VERSION >= 0x060000
	QwtSymbol* symbol = 0;
	if (symbolName == "UserSymbol")
		symbol = new QwtSymbol(QwtSymbol::NoSymbol, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Ellipse")
		symbol = new QwtSymbol(QwtSymbol::Ellipse, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Rect")
		symbol = new QwtSymbol(QwtSymbol::Rect, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Diamond")
		symbol = new QwtSymbol(QwtSymbol::Diamond, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Triangle")
		symbol = new QwtSymbol(QwtSymbol::Triangle, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "LTriangle")
		symbol = new QwtSymbol(QwtSymbol::LTriangle, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "RTriangle")
		symbol = new QwtSymbol(QwtSymbol::RTriangle, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Cross")
		symbol = new QwtSymbol(QwtSymbol::Cross, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "XCross")
		symbol = new QwtSymbol(QwtSymbol::XCross, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "HLine")
		symbol = new QwtSymbol(QwtSymbol::HLine, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "VLine")
		symbol = new QwtSymbol(QwtSymbol::VLine, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Star1")
		symbol = new QwtSymbol(QwtSymbol::Star1, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Star2")
		symbol = new QwtSymbol(QwtSymbol::Star2, Qt::NoBrush, pen, QSize(10, 5));
	else if (symbolName == "Hexagon")
		symbol = new QwtSymbol(QwtSymbol::Hexagon, Qt::NoBrush, pen, QSize(10, 5));
	else
		symbol = new QwtSymbol(QwtSymbol::NoSymbol, Qt::NoBrush, pen, QSize(10, 5));
#else
		QwtSymbol symbol;
		if (symbolName == "UserSymbol")
			symbol = QwtSymbol(QwtSymbol::NoSymbol, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Ellipse")
			symbol = QwtSymbol(QwtSymbol::Ellipse, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Rect")
			symbol = QwtSymbol(QwtSymbol::Rect, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Diamond")
			symbol = QwtSymbol(QwtSymbol::Diamond, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Triangle")
			symbol = QwtSymbol(QwtSymbol::Triangle, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "LTriangle")
			symbol = QwtSymbol(QwtSymbol::LTriangle, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "RTriangle")
			symbol = QwtSymbol(QwtSymbol::RTriangle, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Cross")
			symbol = QwtSymbol(QwtSymbol::Cross, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "XCross")
			symbol = QwtSymbol(QwtSymbol::XCross, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "HLine")
			symbol = QwtSymbol(QwtSymbol::HLine, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "VLine")
			symbol = QwtSymbol(QwtSymbol::VLine, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Star1")
			symbol = QwtSymbol(QwtSymbol::Star1, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Star2")
			symbol = QwtSymbol(QwtSymbol::Star2, Qt::NoBrush, pen, QSize(10, 5));
		else if (symbolName == "Hexagon")
			symbol = QwtSymbol(QwtSymbol::Hexagon, Qt::NoBrush, pen, QSize(10, 5));
		else
			symbol = QwtSymbol(QwtSymbol::NoSymbol, Qt::NoBrush, pen, QSize(10, 5));
#endif
	curve->setSymbol(symbol);
}
//__________________________________________________________________________________________________
void Plot::addCurve(const QString & title, const QString & symbol, const bool ySite, int /*offset*/, const QVector<double> &xData, const QVector<double> &yData)
{
	QwtPlotCurve* curve = createCurve(title, symbol, ySite);

#if QWT_VERSION >= 0x060000
	curve->setSamples(xData, yData);
#else
	curve->setData(xData, yData);
#endif
	curve->attach(this);
	dynamic_cast< get::PlotCurve* >(curve)->initializeLegend();
	arrayOfCurves.append(curve);
	setAutoReplot(true);
	replot();
}
//__________________________________________________________________________________________________
void Plot::addCurve(const QString & title, const QString & symbol, const bool yAxisPos, const int size, const double* xData, const double* yData)
{
	QwtPlotCurve* curve = createCurve(title, symbol, yAxisPos);

#if QWT_VERSION >= 0x060000
	curve->setSamples(xData, yData, size);
#else
	curve->setData(xData, yData, size);
#endif
	curve->attach(this);
	dynamic_cast< get::PlotCurve* >(curve)->initializeLegend();
	arrayOfCurves.append(curve);
	setAutoReplot(true);
	replot();
}
//__________________________________________________________________________________________________
#if QWT_VERSION < 0x060100
void Plot::showCurve(QwtPlotItem *item, bool on)
{
	item->setVisible(on);

	QwtLegendItem *legendItem = qobject_cast<QwtLegendItem *>(legend()->find(item));
	if (legendItem)
		legendItem->setChecked(on);
	replot();
}
#else
void Plot::showCurve(QwtPlotItem * /* item */, bool /* on */)
{
}
#endif
//__________________________________________________________________________________________________
#if QWT_VERSION >= 0x060100
void Plot::showCurve(const QVariant & itemInfo, bool on, int /* index */)
{
	QwtPlotItem* item = infoToItem(itemInfo);
	item->setVisible(on);

	replot();
}
#else
void Plot::showCurve(const QVariant & /* itemInfo */, bool /* on */, int /* index */)
{

}
#endif
//__________________________________________________________________________________________________
int Plot::randNb(int min, int max)
{
	int iSecret = 0;
	/* initialize random seed: */

	/* generate secret number: */
	iSecret = rand() % (max - min) + min;
	return iSecret;
}
//__________________________________________________________________________________________________
void Plot::toggleAllCurves()
{
	// Find out whether at least one curve is visible
	bool isAnyVisible = false;
	QVectorIterator<QwtPlotCurve*> i(arrayOfCurves);
	while (i.hasNext())
	{
		if (i.next()->isVisible())
		{
			isAnyVisible = true;
			break;
		}
	}

	// Show or hide all curves
	i.toFront();
	while (i.hasNext())
	{
		i.next()->setVisible(not isAnyVisible);
	}
}
//__________________________________________________________________________________________________
