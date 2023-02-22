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


CemrgDataInteractor::CemrgDataInteractor(){

}

void CemrgDataInteractor::Initialise(QStringList &options) {
    m_dialog = new QDialog(0,0);
    m_dialog->setWindowTitle("User input");

    m_controls.setupUi(m_dialog);
    m_controls.m_comboBox->addItems(options);

    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(accepted()), m_dialog, SLOT(accept()));
    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(rejected()), m_dialog, SLOT(reject()));
}

CemrgDataInteractor::~CemrgDataInteractor() {
    m_dialog->close();
    m_dialog->deleteLater();
}

void CemrgDataInteractor::AddPoint(mitk::StateMachineAction *, mitk::InteractionEvent *interactionEvent) {
    
    // Show the dialog box
    int dialogCode = m_dialog->exec();
    if (dialogCode == QDialog::Rejected) {
        m_dialog->close();
        m_dialog->deleteLater();
        return; // User clicked cancel, so don't add the point
    }

    QString option = m_controls.m_comboBox->currentText();

    // Add the new point with the label and option
    Superclass::AddPoint(nullptr, interactionEvent);
    mitk::Point3D new_pt = GetLastPoint();
    
    // new_pt->SetStringProperty("option", option.toStdString());
    std::cout << "LABEL, " << option.toStdString() << "(" << new_pt[0] << "," << new_pt[1] << "," << new_pt[2] << ")" << std::endl;
}