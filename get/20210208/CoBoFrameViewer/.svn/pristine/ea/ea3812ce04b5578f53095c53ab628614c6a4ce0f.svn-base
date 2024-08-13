#ifndef get_viewgraph_h_INCLUDED
#define get_viewgraph_h_INCLUDED

#include <QWidget>
#include <QString>

class QwtPlotZoomer;
class QwtPlotPicker;
class QwtPlotPanner;
class QwtSymbol;
class Plot;
#include <qwt_global.h>

class ViewGraph: public QWidget
{
Q_OBJECT
public:
	ViewGraph(QWidget *parent = 0);
	~ViewGraph();
	void addCurve(const int numPoints, const double* xData, const double* yData,
			const QString & title = QString(), const QString & symbolName = QString("UserSymbol"),
			const QColor & color = Qt::black, const bool yAxisLeft = true);
	void clear();
	void validatePlot();
	void setTitle(QString title);
	QString getTitle() const;
	void setAxisTitle(int axisId, QString title);
	void updateAxes();
	int nbCurve();
public slots:
	void toggleAllCurves();
private:
	Plot *d_plot;
signals:
	void statusMessage(const QString &);
private Q_SLOTS:
	void moved(const QPoint &);
	void selected();
#ifndef QT_NO_PRINTER
	void print();
#endif
protected:
	QString buidBaseName() const;
	bool exportDocument(const QString & fileName);
private slots:
	void exportDocument();
	void enableZoomMode(bool);
	void zoomOut();
private:
	void showInfo(QString text = QString::null);
	// Widget panneau affichage de courbes
	QwtPlotZoomer *d_zoomer[2];
	QwtPlotPicker *d_picker;
	QwtPlotPanner *d_panner;
};

#endif // get_viewgraph_h_INCLUDED
