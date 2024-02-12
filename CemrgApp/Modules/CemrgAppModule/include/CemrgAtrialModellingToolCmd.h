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
 * Four Chamber Tools (inherits CemrgCommandLine)
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgAtrialModellingToolCmd_h
#define CemrgAtrialModellingToolCmd_h

#include <QString>
#include <QDir>

#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkCommon.h>

// VTK
#include <vtkSmartPointer.h>
#include <vtkConnectivityFilter.h>
#include <vtkIdList.h>
#include <vtkRegularPolygonSource.h>

// ITK
#include <itkImage.h>
#include <itkImageRegionIterator.h>

#include "FourChamberCommon.h"
#include "CemrgCommandLine.h"

enum AtialToolkitMode
{
    UAC,
    FIBREMAP,
    SCALARMAP,
    LATFIELD,
    STIM,
    LABELS,
    VIS,
    PARFILE
};
enum AtrialToolkitParamfiles
{
    RAA,
    PV2_4Ch,
    PA_4Ch,
    LS,
    single_LR_P,
    single_LR_A,
    EE_RA,
    PV4,
    LAA_4Ch,
    alpha,
    PV1,
    PV3_4Ch,
    single_UD_P,
    alpha_RA,
    PA_4Ch_RA,
    EE,
    LS_4Ch_RA,
    LS_4Ch,
    LAA,
    PV4_4Ch,
    PV3,
    beta,
    beta_RA,
    PA_RA,
    PV2,
    PA,
    PV1_4Ch,
    single_UD_A,
    LS_RA
};

enum AtrialToolkitScalarMap {
    SAN, 
    CT, 
    PM, 
    BB
};
class MITKCEMRGAPPMODULE_EXPORT CemrgAtrialModellingToolCmd : public CemrgCommandLine
{

public:
    mitkClassMacro(CemrgAtrialModellingToolCmd, CemrgCommandLine)
    // [uac|fibremap|scalarmap|latfield|stim|labels|vis|getparfile

    void SetModeOfOperation(AtrialToolkitMode atkm);

    QString GetParfile(QtrialToolkitParamfiles atkm, QString meshname="");
    QString GetParfile(QString parfileName, QString meshname = ""); // meshname="" means use default name

    QString UniversalAtrialCoordinates(QString stage, QString atrium, QString layer, QString fibre, QString meshname, QStringList tags, QStringList landmarks, bool fourch=false, bool noraa=false, int scale=1000);
    QString FibreMapping(QString atrium, QString layer, QString fibre, QString meshname, bool msh_endo_epi, QString output, bool fourch, QString tags, QString biproj);
    QString ScalarMapping(QString atrium, QString meshname, bool bb, AtrialToolkitScalarMap scalatSuffix);

    // Helper functions
    QStringList PrepareDockerExecution(QString &executableName, QString &outAbsolutePath, AtrialToolkitMode atkm);

    // inlines 
    inline void SetDockerTag(QString tag) { dockerTag = tag; };
    inline void SetVolume(QString volume) { dataVolume = volume; home(dataVolume); }
    inline bool IsVolumeSet() { return dataVolume.isEmpty(); };

protected:
private:
    QString dockerTag; 
    QString dataVolume; // directory mounted in /data inside container
    QDir home;
    QString mode; // mode of operation of container
};
#endif // CemrgAtrialModellingToolCmd_h