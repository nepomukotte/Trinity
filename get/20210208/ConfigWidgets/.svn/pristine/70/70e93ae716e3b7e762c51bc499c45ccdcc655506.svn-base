/**
 * @file ConfigSetWidget.cpp
 * @date 29 août 2013
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the ConfigWidgets software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "ConfigSetWidget.h"
#include <CCfg/Control/ConfigEditorController.h>
#include <CCfg/Ui/ConfigSetWidget.ui.h>
#include <CCfg/Model/ConfigSetModel.h>
#include <QFileDialog>
#include <QToolButton>
#include <QModelIndex>
#include <QtGlobal> // for QT_VERSION
#include <utl/Logging.h>

namespace CCfg
{
namespace Ui
{

ConfigSetWidget::ConfigSetWidget(QWidget* parent)
	: Utils::Ui::DesignedWidget<Ui_ConfigSetWidget, QWidget>(parent),
	  model_(0),
	  mapper()
{
	setEnabled(false);
	connect(ui_->lineEditRepository, SIGNAL(textEdited(const QString&)),
			this, SLOT(onRepositoryEdited(const QString&)));
	connect(ui_->buttonSetFileRepo, SIGNAL(clicked()),
			this, SLOT(onSelectFileRepositoryClicked()));
	connect(ui_->treeViewConfigSet, SIGNAL(doubleClicked(const QModelIndex&)),
			this, SLOT(onItemDoubleClicked(const QModelIndex&)));
}

ConfigSetWidget::~ConfigSetWidget()
{
}

void ConfigSetWidget::setModel(Model::ConfigSetModel* model)
{
	model_ = model;
	setEnabled(model_ != 0);
	if(model == 0)
	{
		return; //TODO: maybe we should also clear the widget
	}
	connect(model, SIGNAL(repositoryChanged()), this, SLOT(onRepositoryChanged()));
	connect(model, SIGNAL(configBindingChanged(const QString&, const QString&, const QString&)),
			this,  SLOT(onConfigBindingChanged(const QString&, const QString&, const QString&)));
	onRepositoryChanged();

	ui_->treeViewConfigSet->setModel(model);
	ui_->treeViewConfigSet->expandAll();

	mapper.reset(new QSignalMapper);
	const Model::ConfigSetModel::ConfigSlotList& configs = model->configSlots();
	for(size_t i = 0; i < configs.size(); ++i)
	{
		QToolButton* button = new QToolButton(ui_->treeViewConfigSet);
#if QT_VERSION < 0x050000
		button->setText(QApplication::translate("ConfigSetWidget", "Set...", 0, QApplication::UnicodeUTF8));
#else
		button->setText(QApplication::translate("ConfigSetWidget", "Set...", 0));
#endif

		const Model::ConfigSetModel::ConfigItem& cfgItem = configs[i].config();
		if(!cfgItem.url().isEmpty())
		{
			onConfigBindingChanged(
				QString::fromStdString(cfgItem.slot().id()),
				QString(),
				cfgItem.url());
		}

		connect(button, SIGNAL(clicked()), mapper.get(), SLOT(map()));
		mapper->setMapping(button, i);

		ui_->treeViewConfigSet->setIndexWidget(
				model->index(ButtonColumn, cfgItem),
				button);
	}

	connect(mapper.get(), SIGNAL(mapped(int)), this, SLOT(onTreeButtonClicked(int)));

	QHeaderView* header = ui_->treeViewConfigSet->header();
#if QT_VERSION < 0x050000
	header->setResizeMode(NameColumn, QHeaderView::Stretch);
	header->setResizeMode(ButtonColumn, QHeaderView::ResizeToContents);
#else
	header->setSectionResizeMode(NameColumn, QHeaderView::Stretch);
	header->setSectionResizeMode(ButtonColumn, QHeaderView::ResizeToContents);
#endif
}

void ConfigSetWidget::onRepositoryChanged()
{
	ui_->lineEditRepository->setText(model_->repository());
	ui_->treeViewConfigSet->setEnabled(!model_->repository().isEmpty());
	ui_->treeViewConfigSet->update();
}

void ConfigSetWidget::onConfigBindingChanged(const QString& /*configslotId*/, const QString& /*oldUrl*/, const QString& newUrl)
{
	emit showConfig(newUrl);
}

void ConfigSetWidget::onTreeButtonClicked(int i)
{
	if(i >= 0 and static_cast<size_t>(i) < model_->configSlots().size())
	{
		QString configUrl =
			QFileDialog::getOpenFileName(this,
				tr("Bind %1 to")
					.arg(QString::fromStdString(model_->configSlots()[i].id())),
				model_->repository(),
				CCfg::Control::ConfigEditorController::FILE_EXTENSIONS_FILTER);

		if(!configUrl.isEmpty())
		{
			model_->setUrl(i, configUrl);
		}
	}
}

void ConfigSetWidget::onSelectFileRepositoryClicked()
{
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::Directory);
	dialog.setDirectory(model_->repository());
	dialog.setNameFilter(CCfg::Control::ConfigEditorController::FILE_EXTENSIONS_FILTER);

	if(dialog.exec() and dialog.selectedFiles().size() == 1)
	{
		model_->setRepository(dialog.selectedFiles()[0]);
	}
}

void ConfigSetWidget::onRepositoryEdited(const QString& text)
{
	model_->setRepository(text);
}

void ConfigSetWidget::onItemDoubleClicked(const QModelIndex& index)
{
	Model::ConfigSetModel::ConfigItem* item = model_->indexToItem<Model::ConfigSetModel::ConfigItem>(index);
	if(item == 0)
	{
		item = &model_->item<Model::ConfigSetModel::ConfigSlot>(index).config();
	}
	const QString& url = item->url();
	if(!url.isEmpty())
	{
		emit showConfig(url);
	}
}

} /* namespace Ui */
} /* namespace CCfg */
