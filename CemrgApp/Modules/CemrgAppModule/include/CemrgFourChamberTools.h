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

#ifndef CemrgFourChamberTools_h
#define CemrgFourChamberTools_h

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberTools : public CemrgCommandLine {
    public:

        mitkClassMacro(CemrgFourChamberTools, CemrgCommandLine)
        itkFactorylessNewMacro(Self)
        itkCloneMacro(Self)

        bool CheckCarpDirectory();

        inline void SetCarpDirectory(QString carpDir) { _carp_dir = carpDir; };
        inline void SetDockerImageOpenCarp() { _docker_image = "cemrgapp/opencarp"; };

        inline QString CARP_DIR() { return _carp_dir; };

    protected:

    private: 
        QString _carp_dir = "";

};
#endif // CemrgFourChamberTools_h