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

// Qt
#include <QDir>
#include <QDirIterator>
#include <QStringList>

#include "CemrgFourChamberTools.h"

CemrgFourChamberTools::CemrgFourChamberTools() {

}

bool CemrgFourchamberTools::CheckCarpDirectory() {
    bool success = false;
    if (_carp_dir == "") {
        MITK_ERROR << "CARP directory not set!";
    }

    QDir carpDir(_carp_dir);
    if (!carpDir.exists()) {
        MITK_ERROR << "CARP directory does not exist!";
    }
    
    // iterate carpDir and check for required files
    QStringList requiredFiles = {"mguvc", "GlVTKConvert", "GlRuleFibres", "GlElemCenters", "carp.pt", "igbextract"};
    QDiriterator it(carpDir, QDirIterator::Subdirectories);
    while (it.hasNext()) {
        QString file = it.next();
        if (requiredFiles.contains(file)) {
            requiredFiles.removeOne(file);
        }
    }

    success = requiredFiles.isEmpty(); // means all files were found

    return success;
}