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

#ifndef CemrgFourChamberCmd_h
#define CemrgFourChamberCmd_h

#include <QString>

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

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberCmd : public CemrgCommandLine {
    
    public:

        mitkClassMacro(CemrgFourChamberCmd, CemrgCommandLine)

        static bool CheckCarpDirectory(QString dir);
        inline bool CheckCarpDirectory() { return CheckCarpDirectory(_carp_dir); };
        inline void SetCarpless(bool value) { carpless = value; };

        bool ExecuteCarpless(QString executableName, QStringList arguments, QString outputPath, bool isOutputFile = true);

        bool CalculateUvcs(QString base_dir, FourChamberSubfolders fourch_sdirs, QString mesh_sdir, QString meshname, QString input_tags_parfile, QString etags_sdir, QString apex_sdir);
        QString CalculateEndoToEpiLaplace(QString base_dir, FourChamberSubfolders fourch_sdirs, QString meshname, QString endo_surf, QString epi_surf, QString parfile, QString outdir);

        // CARP Utilities
        bool ExecuteMguvc(QString directory, QString model_name, QString input_model, QString output_model, QString np, QString tags_file, QString output_dir, bool laplace_solution, bool custom_apex, QString id_solver = "");
        bool ExecuteGlVTKConvert(QString directory, QString model, QStringList n_list, QString output_dir, bool trim_names = false);
        bool ExecuteGlRuleFibers(VentricularFibresParams vfib);
        bool ExecuteGlRuleFibers(QString directory, QString m, QString type, QString a, QString e, QString l, QString r, double a_endo, double a_epi, double b_endo, double b_epi, QString output_pre);
        bool ExecuteCarp_Pt(QString directory, QString meshname, QString par_sdir, QString parfile, QStringList stim_files, QString output_dir);
        bool ExecuteIgbextract(QString directory, QString sdir, double small_f, double big_F, QString outname="", QString name="phie.igb");
        bool ExecuteGlElemCenters(QString meshPath, QString outputPath); 

        // CARP binaries getters
        inline void SetCarpDirectory(QString carpDir) { _carp_dir = carpDir; };
        inline QString CARP_DIR(QString subpath) { return _carp_dir + "/" + subpath; };

        // Docker calling cemrg/4ch
        inline void SetDockerImageFourch(QString tag = "latest") { SetDockerImage("cemrg/4ch:" + tag); };
        QString DockerExtractSurfaces(QString baseDirectory, QString parFolder, QString inputTagsFilename, QString apexSeptumFolder, QString meshname);
        QString DockerCorrectFibres(QString baseDirectory, QString meshname);
        QString DockerMeshtoolGeneric(QString directory, QString command, QString subcommand, QStringList arguments, QString expectedOutput);

        inline QString mguvc() { return CARP_DIR("mguvc"); };
        inline QString GlVTKConvert() { return CARP_DIR("GlVTKConvert"); };
        inline QString GlRuleFibers() { return CARP_DIR("GlRuleFibers"); };
        inline QString GlElemCenters() { return CARP_DIR("GlElemCenters"); };
        inline QString carp_pt() { return CARP_DIR("carp.pt"); };
        inline QString igbextract() { return CARP_DIR("igbextract"); };

    protected:
        

    private:
        QString _carp_dir = "";
        bool carpless;
};
#endif // CemrgFourChamberCmd_h
