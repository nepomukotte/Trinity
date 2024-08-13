#include "Console.h"
#include <QTextEdit>
#include <QCloseEvent>
#include <QVBoxLayout>
#include <QCoreApplication>

Console::Console(QWidget * parent) : QWidget(parent,
		Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint)
{
	textEdit = new QTextEdit(this);
	textEdit->setReadOnly(true);
	textEdit->setLineWrapMode(QTextEdit::NoWrap);
	textEdit->document()->setMaximumBlockCount(10000);

	QBoxLayout* layout = new QVBoxLayout;
	layout->addWidget(textEdit);
	setLayout(layout);
}

void Console::clear()
{
	textEdit->clear();
}

void Console::displayMessage(const QString & msg)
{
	textEdit->append(msg);
	QCoreApplication::processEvents();
}

/** Ignores user attempts to close the window. */
void Console::closeEvent(QCloseEvent* event)
{
	event->ignore();
}


