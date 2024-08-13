#ifndef UTILS_UI_FILECHOOSER_H_INCLUDED
#define UTILS_UI_FILECHOOSER_H_INCLUDED

#include <QWidget>
#include <QDir>

// forward declarations
class QLineEdit;
class QPushButton;

namespace Utils
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This widget enables to choose a file by entering a path in an editor or by
 * selecting it in a dialog.
 *
 * @author fchateau
 */
class FileChooser : public QWidget
{
	Q_OBJECT

	Q_ENUMS(SelectionMode)
	Q_ENUMS(PathMode)

	Q_PROPERTY(SelectionMode selectionMode READ getSelectionMode WRITE setSelectionMode)
	Q_PROPERTY(PathMode pathMode READ getPathMode WRITE setPathMode)
	Q_PROPERTY(QString fileName READ getFileName WRITE setFileName)
	Q_PROPERTY(QString filter READ getFilter WRITE setFilter)
	Q_PROPERTY(QDir baseDirectory READ getBaseDirectory WRITE setBaseDirectory)

	Q_PROPERTY(QString text READ getFileName WRITE setFileName)

public:
	// enumerations
	enum SelectionMode {FileToLoad, FileToSave, Directory};
	enum PathMode {Absolute, Relative};

	FileChooser(QWidget* parent = 0);
	virtual ~FileChooser();

	SelectionMode getSelectionMode() const;
	PathMode getPathMode() const;
	QString getFileName() const;
	const QString& getFilter() const;
	const QDir& getBaseDirectory() const;

public slots:
	void setSelectionMode(SelectionMode m);
	void setPathMode(PathMode m);
	void setFileName(const QString& fileName);
	void setFilter(const QString& filter);
	void setBaseDirectory(const QDir& dir);

signals:
	void fileNameChanged(const QString& fileName);

private slots:
	void chooseFile();

private:
	FileChooser(const FileChooser& rvalue);
	FileChooser& operator=(const FileChooser& rvalue);

	QLineEdit* pathEditor;
	QPushButton* button;

	QString filter;	///< filters the authorized file extensions in selection dialogs.
	QDir baseDirectory;	///< the base directory for file dialogs and for computing relative paths.
	SelectionMode selectionMode;	///< indicates the type of file to select.
	PathMode pathMode; ///< indicates whether paths must be kept absolute or converted to relative path.
};

}
}

#endif // UTILS_UI_FILECHOOSER_H_INCLUDED
