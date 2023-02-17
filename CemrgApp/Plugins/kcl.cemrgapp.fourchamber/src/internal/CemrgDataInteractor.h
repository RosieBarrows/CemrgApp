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

#include <mitkPointSetInteractor.h>
#include <mitkPointSet.h>
#include <mitkDataNode.h>

#include "ui_FourChamberViewPointLabelSelect.h"

class CemrgDataInteractor : public mitk::PointSetInteractor {

    Q_OBJECT

public:
    CemrgDataInteractor(const QStringList& options);
    ~CemrgDataInteractor() override;

protected:
    void OnAddPoint(StateMachineAction *, InteractionEvent *interactionEvent) override;

private:
    Ui::FourChamberViewPointLabelSelect m_controls;
    QDialog *m_dialog;
};