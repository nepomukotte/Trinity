#ifndef CCFG_UI_CLOSABLETABWIDGET_H_INCLUDED
#define CCFG_UI_CLOSABLETABWIDGET_H_INCLUDED

#include <QTabWidget>

namespace CCfg
{
namespace Ui
{

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is an extended TabWidget providing a close button.
 * It also improves original the QTabView interface by providing signals for tab
 * insertion and removal, by fixing missing currentChanged signal emission in
 * some cases (removal of a tab replaced by another of the same index).
 *
 * @author fchateau
 */
class ClosableTabWidget : public QTabWidget
{
	Q_OBJECT
public:
	ClosableTabWidget(QWidget* parent = 0);
	virtual ~ClosableTabWidget();

	// @{ Removal related methods.
	bool removeTab(int index, bool interactive = false);
	bool removeTab(QWidget* widget, bool interactive = false);
	bool removeAll(bool interactive = false);
	// @}

	using QTabWidget::currentChanged;
signals:
	void currentChanged(QWidget* widget);
	void tabCreated(int index);
	void tabClosing(int index);
	void tabClosed(int index);
	void canClose(int index, bool* allowed);
	void canCloseAll(bool* allowed);

protected:
	// @{ Override of mutation notification methods
	virtual void tabInserted(int index);
	virtual void tabRemoved(int index);
	// @}

protected slots:
	void onCurrentChanged(int index);
	void onCloseTab();

private:
	ClosableTabWidget(const ClosableTabWidget& right);
	ClosableTabWidget& operator=(const ClosableTabWidget& right);
};

}
}

#endif // CCFG_UI_CLOSABLETABWIDGET_H_INCLUDED
