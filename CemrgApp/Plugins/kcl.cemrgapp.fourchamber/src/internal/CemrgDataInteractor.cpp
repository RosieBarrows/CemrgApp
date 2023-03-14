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
#include "mitkInteractionPositionEvent.h"
#include <mitkPointOperation.h>
//
#include "mitkBaseRenderer.h"
#include "mitkDispatcher.h"
#include <mitkPropertyList.h>

#include <QFileInfo>
#include <QString>

#include "CemrgCommonUtils.h"

CemrgDataInteractor::CemrgDataInteractor(){

}

void CemrgDataInteractor::Initialise(QStringList &options, QString path_to_file) {
    m_dialog = new QDialog(0,0);
    m_dialog->setWindowTitle("User input");

    m_controls.setupUi(m_dialog);
    m_controls.m_comboBox->addItems(options);

    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(accepted()), m_dialog, SLOT(accept()));
    QObject::connect(m_controls.m_ok_cancel_button, SIGNAL(rejected()), m_dialog, SLOT(reject()));

    path_to_json = path_to_file;
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
    m_controls.m_comboBox->removeItem(m_controls.m_comboBox->currentIndex());

    // Add the new point with the label and option
    Superclass::AddPoint(nullptr, interactionEvent);
    mitk::Point3D new_pt = GetLastPoint(interactionEvent);

    QFileInfo fi(path_to_json);
    bool fileExists = fi.exists();

    QString pt_str = QString::number(new_pt[0]) + "," + QString::number(new_pt[1]) + "," + QString::number(new_pt[2]);
    MITK_INFO << ("POINT: " + pt_str).toStdString();

    if (fileExists) {
        MITK_INFO(CemrgCommonUtils::ModifyJSONFile(fi.absolutePath(), fi.fileName(), option, pt_str, "array")) << "Added point to file";
    }
}

mitk::Point3D CemrgDataInteractor::GetLastPoint(mitk::InteractionEvent *interactionEvent) {
    mitk::InteractionPositionEvent *positionEvent = dynamic_cast<mitk::InteractionPositionEvent*>(interactionEvent);
    mitk::Point3D res;
    if(positionEvent == nullptr){
        std::cout << "position event not captured" << std::endl;
        res.Fill(0.0);
    } else {
        res = positionEvent->GetPositionInWorld();
    }
    return res;
}