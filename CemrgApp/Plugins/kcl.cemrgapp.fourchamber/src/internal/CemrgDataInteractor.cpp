/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*=========================================================================
 *
 * Cemrg Data Interactor
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * jose.solislemus@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include "CemrgDataInteractor.h"
#include "mitkInternalEvent.h"
#include "mitkMouseMoveEvent.h"
#include "mitkRenderingManager.h"
#include <mitkPointOperation.h>
//
#include "mitkBaseRenderer.h"
#include "mitkDispatcher.h"
#include <mitkPropertyList.h>

CemrgDataInteractor::CemrgDataInteractor(const QStringList &options) {
    m_dialog = new QDialog;
    m_dialog.setWindowTitle("User input");

    m_controls.setupUi(m_dialog);

    m_comboBox = new QComboBox;
    m_comboBox->addItems(options);

    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(accepted()), m_dialog, SLOT(accept()));
    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(rejected()), m_dialog, SLOT(reject()));
}

CemrgDataInteractor::~CemrgDataInteractor() {
    delete m_dialog;
}

void CemrgDataInteractor::OnAddPoint(StateMachineAction *, InteractionEvent *interactionEvent) override {
    if (interactionEvent->GetSender()->GetNumberOfSelected() > 0){
        return; // Only handle events where exactly one point is selected
    }

    // Show the dialog box
    if (m_dialog->exec() == QDialog::Rejected){
        return; // User clicked cancel, so don't add the point
    }

    QString option = m_comboBox->currentText();

    // Add the new point with the label and option
    Superclass::OnAddPoint(nullptr, interactionEvent);
    auto newPointSet = interactionEvent->GetSender()->GetPointSet();
    auto newPoint = newPointSet->GetPoint(newPointSet->GetSize() - 1);

    newPoint->SetStringProperty("option", option.toStdString());
}