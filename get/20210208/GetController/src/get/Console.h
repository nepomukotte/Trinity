#ifndef get_Console_h_INCLUDED
#define get_Console_h_INCLUDED

#include <QWidget>
class QCloseEvent;
class QTextEdit;

class Console : public QWidget
{
	Q_OBJECT
public:
	Console(QWidget * parent = 0);
public slots:
	void clear();
	void displayMessage(const QString & msg);
protected:
	virtual void closeEvent(QCloseEvent* event);
private:
	QTextEdit* textEdit;
};

#endif // get_Console_h_INCLUDED
