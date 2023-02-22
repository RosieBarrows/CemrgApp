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

#include <mitkPointSet.h>
#include <mitkDataNode.h>
#include <mitkDataInteractor.h>
#include <mitkPointSetDataInteractor.h>
#include <QInputDialog>

#include "ui_FourChamberViewPointLabelSelect.h"

class CemrgDataInteractor : public mitk::PointSetDataInteractor {

public:
    mitkClassMacro(CemrgDataInteractor, mitk::PointSetDataInteractor)
    itkFactorylessNewMacro(Self)
    CemrgDataInteractor();

    void Initialise(QStringList &options);
    inline mitk::Point3D GetLastPoint() { return Superclass::m_LastPoint; };

protected:
    virtual ~CemrgDataInteractor();

    void AddPoint(mitk::StateMachineAction *, mitk::InteractionEvent *interactionEvent) override;

private:
    Ui::FourChamberViewPointLabelSelect m_controls;
    QDialog *m_dialog;
};
