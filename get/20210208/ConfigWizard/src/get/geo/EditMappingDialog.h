/**
 * @file EditMappingDialog.h
 * @date Jan 28, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: EditMappingDialog.h 1104 2014-01-28 10:07:20Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the ConfigWizard software project.
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

#ifndef get_geo_EditMappingDialog_h_INCLUDED
#define get_geo_EditMappingDialog_h_INCLUDED

#include "get/geo/EditMappingDialog.ui.h"
#include <Utils/Ui/DesignedWidget.hpp>
#include <QDialog>

namespace get
{
namespace geo
{
class GeometryModel;
//__________________________________________________________________________________________________
class EditMappingDialog : public Utils::Ui::DesignedWidget< Ui_EditMappingDialog, QDialog >
{
	Q_OBJECT
public:
	EditMappingDialog(QWidget* parentWidget = 0);
	QSize sizeHint() const;
	void setGeometryModel(GeometryModel* model);
	GeometryModel* geometryModel() const;
public slots:
	void on_addButton_clicked(bool);
	void on_removeButton_clicked(bool);
	void on_clearButton_clicked(bool);
	void on_importButton_clicked(bool);
	void on_exportButton_clicked(bool);
	void accept();
};
//__________________________________________________________________________________________________
} /* namespace geo */
} /* namespace get */
#endif /* get_geo_EditMappingDialog_h_INCLUDED */
