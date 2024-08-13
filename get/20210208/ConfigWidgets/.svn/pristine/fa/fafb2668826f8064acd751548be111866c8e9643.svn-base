#include "FileChooser.h"
#include <QLineEdit>
#include <QFileDialog>
#include <QPushButton>
#include <QHBoxLayout>

namespace Utils
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// Utils::Ui::FileChooser class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
FileChooser::FileChooser(QWidget* parent)
	: QWidget(parent),
	pathEditor(new QLineEdit(this)),
	button(new QPushButton(this)),
	filter(),
	baseDirectory(QDir::currentPath()),
	selectionMode(FileToLoad),
	pathMode(Absolute)
{
	button->setText("...");
	button->setFixedWidth(button->fontMetrics().width(" ... "));

	QHBoxLayout *layout = new QHBoxLayout(this);
    layout->setMargin(0);
    layout->addWidget(pathEditor);
    layout->addWidget(button);
	setFocusProxy(pathEditor);

	connect(pathEditor, SIGNAL(textChanged(QString)), this, SIGNAL(fileNameChanged(QString)));
	connect(button, SIGNAL(clicked()), this, SLOT(chooseFile()));
}

/**
 * Destructor.
 */
FileChooser::~FileChooser()
{
}

/**
 * Returns the selection mode.
 */
FileChooser::SelectionMode FileChooser::getSelectionMode() const
{
	return selectionMode;
}

/**
 * Returns the path mode.
 */
FileChooser::PathMode FileChooser::getPathMode() const
{
	return pathMode;
}

/**
 * Returns the file name in the editor.
 */
QString FileChooser::getFileName() const
{
	return pathEditor->text();
}

/**
 * Returns the filter string for file selection, which allows to filter file extensions.
 */
const QString& FileChooser::getFilter() const
{
	return filter;
}

/**
 * Defines the base directory for file dialogs and for computing relative paths.
 */
const QDir& FileChooser::getBaseDirectory() const
{
	return baseDirectory;
}

/**
 * Defines the selection mode (file to load, file to save, or directory).
 * The default mode is FileToLoad.
 */
void FileChooser::setSelectionMode(SelectionMode m)
{
	selectionMode = m;
}

/**
 * Defines the path mode (whether paths are converted to relative or absolute paths).
 */
void FileChooser::setPathMode(PathMode m)
{
	pathMode = m;
}

/**
 * Defines the file name in the editor.
 */
void FileChooser::setFileName(const QString& fileName)
{
	if(pathMode == Absolute)
	{
		pathEditor->setText(baseDirectory.absoluteFilePath(fileName));
	}
	else
	{
		pathEditor->setText(baseDirectory.relativeFilePath(fileName));
	}
}

/**
 * Defines the filter string to filter some extensions in the file selection dialog.
 * It must use the same conventions as the QFileDialog.
 */
void FileChooser::setFilter(const QString& filter)
{
	this->filter = filter;
}

/**
 * Defines the base directory for file dialogs and for computing relative paths.
 */
void FileChooser::setBaseDirectory(const QDir& dir)
{
	baseDirectory = dir;
}

/**
 * Display a dialog to choose a file (or directory, according to the mode).
 */
void FileChooser::chooseFile()
{
    QString fileName;
    QDir initFileName(baseDirectory.absoluteFilePath(getFileName()));
    switch(selectionMode)
    {
    case FileToLoad:
    	fileName = QFileDialog::getOpenFileName(this,
    			QString(),
    			initFileName.absolutePath(),
    			filter);
    	break;
    case FileToSave:
    	fileName = QFileDialog::getSaveFileName(this,
    			QString(),
    			initFileName.absolutePath(),
    			filter);
    	break;
    case Directory:
    	fileName = QFileDialog::getExistingDirectory(this,
    			QString(),
    			initFileName.absolutePath());
    	break;
    }

    if(not fileName.isEmpty())
    {
    	setFileName(fileName);
    }
    pathEditor->setFocus();
}

}
}

