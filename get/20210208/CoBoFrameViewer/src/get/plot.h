#ifndef _PLOT_H_
#define _PLOT_H_ 

#include <qwt_plot.h>
#include <qwt_symbol.h>
#include <QVector>

class QwtPlotCurve;
class QwtPlotMarker;
class QWtPlotItem;

class Plot: public QwtPlot
{
    Q_OBJECT
public:
    Plot(QString title,QWidget *parent);
    ~Plot();
    void clear();
public Q_SLOTS:
	void addCurve(const int numPoints, const double* xData, const double* yData,
			const QString & title = QString(), const QString & symbolName = QString("UserSymbol"),
			const QColor & color = Qt::black, const bool yAxisLeft = true);
    void toggleAllCurves();
    int size();
private Q_SLOTS:
	void showCurve(QwtPlotItem *, bool on);
	void showCurve(const QVariant &, bool, int);
private:
    QwtPlotCurve* createCurve(const QString & title, const QString & symbolName, const QColor & color = Qt::black, const bool yAxisLeft = true);
    static QwtSymbol::Style symbolNameToStyle(const QString & symbolName);
    void setCurveSymbol(QwtPlotCurve* curve, const QString & symbolName, const QPen & pen);

    QVector<QwtPlotCurve*> arrayOfCurves;
    QwtPlotMarker *d_marker1;
    QwtPlotMarker *d_marker2;
};

#endif
