#ifndef CCFG_UI_MYCOMBOBOX_H_INCLUDED
#define CCFG_UI_MYCOMBOBOX_H_INCLUDED

#include <QComboBox>
#include <QString>

namespace CCfg
{
namespace Ui
{
/**
 * This is a combo box that has an additionnal "text" property that operates
 * directly on item text rather than item indexes.
 * It is necessary for the ValueEditorCreator, because in a QItemDelegate, the
 * value is set automatically set using the "text" property, but QComboBox does
 * not have this property.
 *
 * @author fchateau.
 */
class MyComboBox : public QComboBox
{
	Q_OBJECT

	Q_PROPERTY(QString text
		READ getText
		WRITE setText
		USER true)
public:
	MyComboBox(QWidget* parent);
	virtual ~MyComboBox();

	QString getText() const;
	void setText(const QString& text);
};

}
}


#endif // CCFG_UI_MYCOMBOBOX_H_INCLUDED
