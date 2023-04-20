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
#include <QString>
#include <QStringList>
#include <QFile>
#include <QFileInfo>

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

bool CemrgFourChamberTools::CalculateUvcs(QString base_dir, FourChamberSubfolders fourch_sdirs, QString mesh_sdir, QString meshname, QString input_tags_parfile, QString etags_sdir, QString apex_sdir) {
    QString etags_folder = base_dir + "/" + fourch_sdirs.PAR + "/" + etags_sdir;
    QString mesh = base_dir + "/" + fourch_sdirs.MESH + "/" + mesh_sdir + "/" + meshname;
    QString surf_folder = base_dir + "/" + fourch_sdirs.UVC; 
    QString surf_folder_la = base_dir + "/" + fourch_sdirs.UVC_LA;
    QString surf_folder_ra = base_dir + "/" + fourch_sdirs.UVC_RA;

    QStringList n_list;
    // constant arguments for mguvc
    const bool LAPSOL = true;
    const bool C_APEX = true;
    const QString ID_SOL = "$MESH/UVC_ek";

    // constant arguments for GlVTKConvert
    const bool TRIM_N = true;

    MITK_INFO << "Calculating UVCs for the BiV mesh...";
    QFile::copy(etags_folder + "/etags.sh", surf_folder + "/BiV/etags.sh");

    // mguvc
    bool mguvc_success = ExecuteMguvc(surf_folder, "BiV/BiV", "biv", "biv", "20", "BiV/etags.sh", "BiV/uvc/", LAPSOL, false);
    
    if (!mguvc_success) {
        MITK_ERROR << "mguvc failed!";
        return false;
    }

    n_list << "BiV/uvc/BiV.uvc_phi.dat";
    n_list << "BiV/uvc/BiV.uvc_z.dat";   
    n_list << "BiV/uvc/BiV.uvc_ven.dat";
    n_list << "BiV/uvc/BiV.uvc_rho.dat";
    bool glvtkconvert_success = ExecuteGlVTKConvert(surf_folder, "BiV/BiV", n_list, "BiV/uvc/uvc", TRIM_N);

    if (!glvtkconvert_success) {
        MITK_ERROR << "GlVTKConvert (BiV/uvc/uvc) failed!";
        return false;
    }

    n_list.clear();

    // GlVTKConvert 
    n_list << "BiV/uvc/BiV.sol_apba_lap.dat";
    n_list << "BiV/uvc/BiV.sol_rvendo_lap.dat";
    n_list << "BiV/uvc/BiV.sol_endoepi_lap.dat";
    n_list << "BiV/uvc/BiV.sol_lvendo_lap.dat";
    glvtkconvert_success = ExecuteGlVTKConvert(surf_folder, "BiV/BiV", n_list, "BiV/uvc/laplace", TRIM_N);

    if (!glvtkconvert_success) {
        MITK_ERROR << "GlVTKConvert (BiV/uvc/laplace) failed!";
        return false;
    }

    n_list.clear();

    MITK_INFO << "Calculating UVCs for the LA mesh...";
    QFile::copy(etags_folder + "/etags_la.sh", surf_folder_la + "/la/etags.sh");

    // mguvc
    mguvc_success = ExecuteMguvc(surf_folder_la, "la/la", "lv", "lv", "20", "la/etags.sh", "la/uvc/", LAPSOL, C_APEX, ID_SOL);

    if (!mguvc_success) {
        MITK_ERROR << "mguvc (1) failed!";
        return false;
    }

    // GlVTKConvert
    n_list <<"la/uvc/la.uvc_phi.dat";
    n_list <<"la/uvc/la.uvc_z.dat";
    n_list <<"la/uvc/la.uvc_ven.dat";
    n_list <<"la/uvc/la.uvc_rho.dat";
    glvtkconvert_success = ExecuteGlVTKConvert(surf_folder_la, "la/la", arguments, "la/uvc/uvc", TRIM_N);

    if (!glvtkconvert_success) {
        MITK_ERROR << "GlVTKConvert (3) failed!";
        return false;
    }

    n_list.clear();

    MITK_INFO << "Calculating UVCs for the RA mesh...";
    QFile::copy(etags_folder + "/etags_ra.sh", surf_folder_ra + "/ra/etags.sh");

    // mguvc
    mguvc_success = ExecuteMguvc(surf_folder_ra + "ra/ra", "lv", "lv", "20", "ra/etags.sh", "ra/uvc/", LAPSOL, C_APEX, ID_SOL);

    if (!mguvc_success) {
        MITK_ERROR << "mguvc (ra/uvc/) failed!";
        return false;
    }

    // GlVTKConvert
    n_list <<"ra/uvc/ra.uvc_phi.dat";
    n_list <<"ra/uvc/ra.uvc_z.dat";
    n_list <<"ra/uvc/ra.uvc_ven.dat";
    n_list <<"ra/uvc/ra.uvc_rho.dat";
    glvtkconvert_success = ExecuteGlVTKConvert(surf_folder_ra, "ra/ra", n_list, "ra/uvc/uvc", TRIM_N);

    MITK_ERROR(!glvtkconvert_success) << "GlVTKConvert (ra/uvc/uvc) failed!";
    
    return glvtkconvert_success;
}

bool CemrgFourChamberTools::ExecuteMguvc(QString directory, QString model_name, QString input_model, QString output_model, QString np, QString tags_file, QString output_dir, bool laplace_solution, bool custom_apex, QString id_solve) {
    if (!id_solve.isEmpty()) {
        arguments << "--ID="+id_solve;
    }

    QString output_path = directory + "/" + output_dir;

    arguments << "--model-name" << directory + "/" + model_name;
    arguments << "--input-model" << input_model << "--output-model" << output_model;
    arguments << "--np" << np;
    arguments << "--tags-file" << directory + "/" + tags_file;
    arguments << "--output-dir" << output_path;
    if (laplace_solution) {
        arguments << "--laplace-solution";
    }
    if (custom_apex) {
        arguments << "--custom-apex";
    }
    
    return this->ExecuteCommand(mguvc(), arguments, output_path, false);
}

bool CemrgFourChamber::ExecuteGlVTKConvert(QString directory, QString model, QStringList n_list, QString output_dir, bool trim_names = false) {
    QString output_path = directory + "/" + output_dir;
    QStringList arguments;

    arguments << "-m" << directory + "/" + model;
    for (int ix = 0; ix < n_list.size(); ix++) {
        arguments << "-n" << directory + "/" + n_list.at(ix);
    }
        
    arguments << "-o" << output_path;
    if (trim_names) {
            arguments << "--trim-names";
    }

    return this->ExecuteCommand(glvtkconvert(), arguments, output_path, false);
}