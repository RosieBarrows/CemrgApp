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

// MITK
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkMorphologicalOperations.h>

// ITK
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

#include "CemrgCommonUtils.h"
#include "CemrgFourChamberCmd.h"

bool CemrgFourChamberCmd::CheckCarpDirectory(QString dir) {
    bool success = false;
    if (dir == "") {
        MITK_ERROR << "CARP directory not set!";
        return success;
    }

    QDir carpDir(dir);
    if (!carpDir.exists()) {
        MITK_ERROR << "CARP directory does not exist!";
        return success;
    }
    
    // iterate carpDir and check for required files
    QStringList requiredFiles = {"mguvc", "GlVTKConvert", "GlRuleFibres", "GlElemCenters", "carp.pt", "igbextract"};
    QDirIterator it(carpDir, QDirIterator::Subdirectories);
    while (it.hasNext()) {
        QString file = it.next();
        for (int ix = 0; ix < requiredFiles.size(); ix++) {
            if (file.contains(requiredFiles.at(ix))) {
                requiredFiles.removeOne(requiredFiles.at(ix));
            }
        }
    }

    if (!requiredFiles.isEmpty()) {
        MITK_WARN << "CARP directory does not contain all required files! Below a list of missing files:";
        for (int ix = 0; ix < requiredFiles.size(); ix++) {
            MITK_WARN << requiredFiles.at(ix).toStdString();
        }
    }
    success = requiredFiles.isEmpty(); // means all files were found

    return success;
}

bool CemrgFourChamberCmd::ExecuteCarpless(QString executableName, QStringList arguments, QString outputPath, bool isOutputFile) {
    bool res = false;
    if (carpless) {
        MITK_WARN << "CARPless mode!";
        MITK_INFO << this->PrintFullCommand(executableName, arguments);
    } else {
        res = this->ExecuteCommand(executableName, arguments, outputPath, isOutputFile);
    }
    return res;
}

bool CemrgFourChamberCmd::CalculateUvcs(QString base_dir, FourChamberSubfolders fourch_sdirs, QString mesh_sdir, QString meshname, QString input_tags_parfile, QString etags_sdir, QString apex_sdir) {
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
    glvtkconvert_success = ExecuteGlVTKConvert(surf_folder_la, "la/la", n_list, "la/uvc/uvc", TRIM_N);

    if (!glvtkconvert_success) {
        MITK_ERROR << "GlVTKConvert (3) failed!";
        return false;
    }

    n_list.clear();

    MITK_INFO << "Calculating UVCs for the RA mesh...";
    QFile::copy(etags_folder + "/etags_ra.sh", surf_folder_ra + "/ra/etags.sh");

    // mguvc
    mguvc_success = ExecuteMguvc(surf_folder_ra, "ra/ra", "lv", "lv", "20", "ra/etags.sh", "ra/uvc/", LAPSOL, C_APEX, ID_SOL);

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

bool CemrgFourChamberCmd::ExecuteMguvc(QString directory, QString model_name, QString input_model, QString output_model, QString np, QString tags_file, QString output_dir, bool laplace_solution, bool custom_apex, QString id_solve) {
    QString output_path = directory + "/" + output_dir;
    QStringList arguments;

    if (!id_solve.isEmpty()) {
        arguments << "--ID="+id_solve;
    }

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
    
    return ExecuteCarpless(mguvc(), arguments, output_path, false);
}

bool CemrgFourChamberCmd::ExecuteGlVTKConvert(QString directory, QString model, QStringList n_list, QString output_dir, bool trim_names) {
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

    return ExecuteCarpless(GlVTKConvert(), arguments, output_path, false);
}

bool CemrgFourChamberCmd::ExecuteCarp_Pt(QString directory, QString meshname, QString par_sdir, QString parfile, QStringList stim_files, QString output_dir) {
    QString output_path = directory + "/" + output_dir;
    QString param_path = directory + "/" + par_sdir + "/" + parfile;

    QStringList arguments;
    arguments << "+F" << param_path;
    arguments << "-simID" << output_path;
    arguments << "-meshname" << directory + "/" + meshname;
    arguments << "-stimulus[0].vtx_file" << stim_files.at(0);
    arguments << "-stimulus[1].vtx_file" << stim_files.at(1);

    return ExecuteCarpless(carp_pt(), arguments, output_path, false);
}

bool CemrgFourChamberCmd::ExecuteIgbextract(QString directory, QString sdir, double small_f, double big_F, QString outname, QString name) {
    QString io_dir = directory + "/" + sdir;
    QString input_igb = io_dir + "/" + name;
    QFileInfo fi(input_igb); 
    
    if (outname.isEmpty()) {
        outname = fi.baseName();
    }

    if (!fi.exists()) {
        MITK_ERROR << ("IGB file [" + input_igb + "] not found").toStdString();
        return false;
    }

    outname += (outname.contains(".")) ? "" : ".dat";
    QString output_path = io_dir + "/" + outname;

    QStringList arguments;
    arguments << input_igb; 
    arguments << "-o" << "ascii";
    arguments << "-f" << QString::number(small_f);
    arguments << "-F" << QString::number(big_F);
    arguments << "-O" << output_path ;

    return ExecuteCarpless(igbextract(), arguments, output_path, true);

}

bool CemrgFourChamberCmd::ExecuteGlElemCenters(QString meshPath, QString outputPath) {
    QStringList arguments;

    arguments << "--meshname" << meshPath;
    arguments << "--output" << outputPath;

    return ExecuteCarpless(GlElemCenters(), arguments, outputPath, true);

}

bool CemrgFourChamberCmd::ExecuteGlRuleFibres(VentricularFibresParams vfib) {
    QStringList arguments;

    arguments << "--meshname" << vfib.meshname();
    arguments << "--type" << vfib.type();
    arguments << "--apex_to_base" << vfib.apex_to_base() ;
    arguments << "--epi" << vfib.epi() ;
    arguments << "--lv" << vfib.lv() ;
    arguments << "--rv" << vfib.rv();
    arguments << "--alpha_endo" << vfib.a_endo();
    arguments << "--alpha_epi" << vfib.a_epi();
    arguments << "--beta_endo" << vfib.b_endo();
    arguments << "--beta_epi" << vfib.b_epi();
    arguments << "--output" << vfib.output();

    return ExecuteCarpless(GlRuleFibres(), arguments, vfib.output(), true);
    
}

bool CemrgFourChamberCmd::ExecuteGlRuleFibres(QString directory, QString m, QString type, QString a, QString e, QString l, QString r, double a_endo, double a_epi, double b_endo, double b_epi, QString output_pre) {
    VentricularFibresParams vfib;
    vfib.SetDirectory(directory); 
    vfib.SetMeshname(m);
    vfib.SetType(type);
    vfib.SetApexToBase(a);
    vfib.SetEpi(e);
    vfib.SetLV(l);
    vfib.SetRV(r);
    vfib.SetAlphaEndo(a_endo);
    vfib.SetAlphaEpi(a_epi);
    vfib.SetBetaEndo(b_endo);
    vfib.SetBetaEpi(b_epi);
    vfib.SetOutput(output_pre);

    return ExecuteGlRuleFibres(vfib);
}