#include "ClosableTabWidget.h"
#include <QToolButton>
#include <QAction>

namespace CCfg
{
namespace Ui
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ClosableTabWidget class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ClosableTabWidget::ClosableTabWidget(QWidget* parent)
	: QTabWidget(parent)
{
	QToolButton* closeButton = new QToolButton(this);
	QAction* closeAction = new QAction(tr("Close tab"), closeButton);
	closeAction->setShortcut(QKeySequence::Close);
	closeAction->setShortcutContext(Qt::WindowShortcut);
	closeAction->setIcon(QIcon(":/tabwidget/closetab.png"));
	closeAction->setToolTip(tr("Close the current tab"));

    closeButton->setDefaultAction(closeAction);
    closeButton->setAutoRaise(true);
    closeButton->adjustSize();

    connect(closeAction, SIGNAL(triggered()), this, SLOT(onCloseTab()));
	connect(this, SIGNAL(currentChanged(int)), this, SLOT(onCurrentChanged(int)));

    setCornerWidget(closeButton);
}

/**
 * Destructor.
 */
ClosableTabWidget::~ClosableTabWidget()
{
}

/**
 * Removes the specified tab.
 * Overloads (but do not overrides) the default method to fix issues with
 * currentChanged(int) signal emission.
 */
bool ClosableTabWidget::removeTab(int index, bool interactive)
{
	bool allowed = true;
	if(interactive)
	{
		emit canClose(index, &allowed);
	}
	if(allowed)
	{
		emit tabClosing(index);
		if(count() > 1)
		{
			if(index == count() - 1)
			{
				setCurrentIndex(index - 1);
			}
			else
			{
				setCurrentIndex(index + 1);
			}
		}
		QTabWidget::removeTab(index);
	}
	return allowed;
}

/**
 * Removes the tab containing the specified widget.
 */
bool ClosableTabWidget::removeTab(QWidget* widget, bool interactive)
{
	int index = indexOf(widget);
	if(index != -1)
	{
		return removeTab(index, interactive);
	}
	else
	{
		return false;
	}
}

/**
 * Removes all tabs.
 */
bool ClosableTabWidget::removeAll(bool interactive)
{
	bool allowed = true;
	if(interactive)
	{
		emit canCloseAll(&allowed);
	}
	if(allowed)
	{

		for(int i = count()-1; i >= 0; --i)
		{
			removeTab(i);
		}
	}
	return allowed;
}

/**
 * Intercepts the signal currentChanged(int) to implement the
 * currentChanged(QWidget*) signal.
 */
void ClosableTabWidget::onCurrentChanged(int index)
{
	QWidget* nextCurrentWidget = widget(index);
	if(nextCurrentWidget != currentWidget())
	{
		emit currentChanged(nextCurrentWidget);
	}
}

/**
 * Reacts to the user action closing a tab.
 */
void ClosableTabWidget::onCloseTab()
{
	removeTab(currentIndex(), true);
}

/**
 * Overrides the base class method to emit signals correctly.
 */
void ClosableTabWidget::tabInserted(int index)
{
	QTabWidget::tabInserted(index);
	cornerWidget()->setEnabled(count() != 0);
	emit tabCreated(index);
	if(count() == 1)
	{
		emit currentChanged(index);
	}
}

/**
 * Overrides the base class method to emit signals.
 */
void ClosableTabWidget::tabRemoved(int index)
{
	QTabWidget::tabRemoved(index);
	cornerWidget()->setEnabled(count() != 0);
	emit tabClosed(index);
}

}
}
