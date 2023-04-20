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

#include <QString>

#include "FourChamberCommon.h"

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberTools : public CemrgCommandLine {
    public:

        mitkClassMacro(CemrgFourChamberTools, CemrgCommandLine)
        itkFactorylessNewMacro(Self)
        itkCloneMacro(Self)

        bool CheckCarpDirectory();

        QString CalculateUvcs(QString base_dir, FourChamberSubfolders fourch_sdirs, QString mesh_sdir, QString meshname, QString input_tags_parfile, QString etags_sdir, QString apex_sdir);

        bool ExecuteMguvc(QString directory, QString model_name, QString input_model, QString output_model, QString np, QString tags_file, QString output_dir, bool laplace_solution, bool custom_apex, QString id_solver = "");
        bool ExecuteGlVTKConvert(QString directory, QString model, QStringList n_list, QString output_dir, bool trim_names = false);

        inline void SetCarpDirectory(QString carpDir) { _carp_dir = carpDir; };
        inline void SetDockerImageOpenCarp() { _docker_image = "cemrgapp/opencarp"; };

        inline QString CARP_DIR() { return _carp_dir; };

        inline QString mguvc() { return _carp_dir + "/mguvc"; };
        inline QString GlVTKConvert() { return _carp_dir + "/GlVTKConvert"; };
        inline QString GlRuleFibres() { return _carp_dir + "/GlRuleFibres"; };
        inline QString GlElemCenters() { return _carp_dir + "/GlElemCenters"; };
        inline QString carp_pt() { return _carp_dir + "/carp.pt"; };
        inline QString igbextract() { return _carp_dir + "/igbextract"; };

    protected:

    private: 
        QString _carp_dir = "";

};
#endif // CemrgFourChamberTools_h