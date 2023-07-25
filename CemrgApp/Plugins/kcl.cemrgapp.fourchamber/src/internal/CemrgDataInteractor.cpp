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
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonValue>
#include <QMessageBox>

#include "CemrgCommonUtils.h"

CemrgDataInteractor::CemrgDataInteractor(){
    m_dialog = nullptr;
    path_to_json = "";
}

void CemrgDataInteractor::Initialise(QStringList &options, QString path_to_file, QString title) {
    m_dialog = std::make_unique<QDialog>();
    m_dialog->setWindowTitle("Select " + title + " Points");

    m_ptset_controls.setupUi(m_dialog.get());
    m_ptset_controls.m_comboBox->addItems(options);

    QObject::connect(m_ptset_controls.m_ok_cancel_button, SIGNAL(accepted()), m_dialog.get(), SLOT(accept()));
    QObject::connect(m_ptset_controls.m_ok_cancel_button, SIGNAL(rejected()), m_dialog.get(), SLOT(reject()));

    path_to_json = path_to_file;
}

CemrgDataInteractor::~CemrgDataInteractor() {
    QObject::disconnect(m_ptset_controls.m_ok_cancel_button, SIGNAL(accepted()), m_dialog.get(), SLOT(accept()));
    QObject::disconnect(m_ptset_controls.m_ok_cancel_button, SIGNAL(rejected()), m_dialog.get(), SLOT(reject()));
    
    m_dialog->close();
    m_dialog->deleteLater();
    m_dialog = nullptr;
}

void CemrgDataInteractor::AddPoint(mitk::StateMachineAction *, mitk::InteractionEvent *interactionEvent) {
    
    // Show the dialog box
    int dialogCode = m_dialog->exec();
    if (dialogCode == QDialog::Rejected) {
        this->~CemrgDataInteractor();
        return; // User clicked cancel, so don't add the point
    }

    if (m_ptset_controls.m_comboBox->count() == 0) {
        QMessageBox::warning(nullptr, "Warning", "No more points to add");
        this->~CemrgDataInteractor();
        return;
    }
    QString option = m_ptset_controls.m_comboBox->currentText();
    m_ptset_controls.m_comboBox->removeItem(m_ptset_controls.m_comboBox->currentIndex());

    // Add the new point with the label and option
    Superclass::AddPoint(nullptr, interactionEvent);
    mitk::Point3D new_pt = GetLastPoint(interactionEvent);

    QFileInfo fi(path_to_json);
    bool fileExists = fi.exists();

    QString pt_str = QString::number(new_pt[0]) + "," + QString::number(new_pt[1]) + "," + QString::number(new_pt[2]);
    MITK_INFO << ("POINT: [" + option + "]" + pt_str).toStdString();

    if (fileExists) {
        MITK_INFO(CemrgCommonUtils::ModifyJSONFile(fi.absolutePath(), fi.fileName(), option, pt_str, "array")) << "Added point to file";
        CheckPoints();
    } else {
        MITK_WARN << "JSON File does not exist";
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

void CemrgDataInteractor::CheckPoints(){
    QFileInfo fi(path_to_json);
    if (!fi.exists()){
        return;
    }

    QJsonObject points = CemrgCommonUtils::ReadJSONFile(fi.absolutePath(), fi.fileName());
    QStringList cylinders = PointIds.CYLINDERS();
    QStringList slicers = PointIds.SLICERS();
    QStringList valve_plains = PointIds.VALVE_PLAINS();

    std::string output = "Selected Points\n\n";

    output += PrintPoints(points, cylinders, "Cylinders") + '\n';
    output += PrintPoints(points, slicers, "Slicers") + '\n';
    output += PrintPoints(points, valve_plains, "Valve Plains") + '\n';

    std::cout << output.c_str();
    QMessageBox::information(nullptr, "Selected Points", output.c_str());
}

std::string CemrgDataInteractor::PrintPoints(QJsonObject json, QStringList keysList, QString title) {
    std::string output = title.toStdString() + "\n";
    for (int ix = 0; ix < keysList.size(); ix++) {
        double arr[3] = {0, 0, 0};
        QString key = keysList.at(ix);
        if (json[key].isUndefined())
        {
            output += key.toStdString() + ": Undefined\n";
            continue;
        }
        output += key.toStdString() + ": ";
        for (int jx = 0; jx < 3; jx++)
        {
            arr[jx] = json[key].toArray().at(jx).toDouble();
        }

        if (arr[0] == 0 && arr[1] == 0 && arr[2] == 0)
        {
            output += "NOT SET\n";
            continue;
        }

        for (int jx = 0; jx < 3; jx++)
        {
            output += QString::number(arr[jx]).toStdString();
            output += (jx < 2) ? ", " : ")\n";
        }
    }

    return output;
}