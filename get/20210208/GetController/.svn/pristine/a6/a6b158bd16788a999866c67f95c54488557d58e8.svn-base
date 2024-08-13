#ifndef _PLOT_H_
#define _PLOT_H_ 

#include <qwt_plot.h>
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
	void addCurve(const QString & title, const QString & symb, const bool ySite, const int offset, const QVector<double> &xData, const QVector<double> & yData);
	void addCurve(const QString & title, const QString & symb, const bool ySite, const int size, const double* xData, const double* yData);
    void toggleAllCurves();
    int size();
private Q_SLOTS:
	void showCurve(QwtPlotItem *, bool on);
	void showCurve(const QVariant &, bool, int);
private:
    QwtPlotCurve* createCurve(const QString & title, const QString & symb, const bool ySite);
    void setCurveSymbol(QwtPlotCurve* curve, const QString & symbolName, const QPen & pen);
    int randNb(int min,int max);

    QVector<QwtPlotCurve*> arrayOfCurves;
    QwtPlotMarker *d_marker1;
    QwtPlotMarker *d_marker2;
};

#endif
