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

#include <mitkIOUtil.h>
#include <mitkImage.h>

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

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberTools : public CemrgCommandLine {
    typedef itk::Image<uint8_t, 3> ImageType;
    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    
    public:

        mitkClassMacro(CemrgFourChamberTools, CemrgCommandLine)
        itkFactorylessNewMacro(Self)
        itkCloneMacro(Self)

        bool CheckCarpDirectory();

        // Segmentation Utilities
        void ExploreLabelsToSplit(mitk::Image::Pointer seg, std::vector<int>& labels); 
        mitk::Image::Pointer SplitLabelsOnRepeat(mitk::Image::Pointer seg, int label);

        QString CalculateUvcs(QString base_dir, FourChamberSubfolders fourch_sdirs, QString mesh_sdir, QString meshname, QString input_tags_parfile, QString etags_sdir, QString apex_sdir);
        QString CalculateEndoToEpiLaplace(QString base_dir, FourChamberSubfolders fourch_sdirs, QString meshname, QString endo_surf, QString epi_surf, QString parfile, QString outdir);

        // CARP Utilities
        bool ExecuteMguvc(QString directory, QString model_name, QString input_model, QString output_model, QString np, QString tags_file, QString output_dir, bool laplace_solution, bool custom_apex, QString id_solver = "");
        bool ExecuteGlVTKConvert(QString directory, QString model, QStringList n_list, QString output_dir, bool trim_names = false);
        bool ExecuteGlRuleFibres(QString directory, VFibresParams vfib, QString output_pre);
        bool ExecuteGlRuleFibres(QString directory, QString m, QString type, QString a, QString e, QString l, QString r, double a_endo, double a_epi, double b_endo, double b_epi, QString output_pre);
        bool ExecuteCarp_Pt(QString directory, QString meshname, QString par_sdir, QString parfile, QStringList stim_files, QString output_dir);
        bool ExecuteIgbextract(QString directory, QString sdir, double small_f, double big_F, QString outname="", QString name="phie.igb");

        // CARP binaries getters
        inline void SetCarpDirectory(QString carpDir) { _carp_dir = carpDir; };
        inline QString CARP_DIR(QString subpath) { return _carp_dir + "/" + subpath; };

        inline QString mguvc() { return CARP_DIR("mguvc"); };
        inline QString GlVTKConvert() { return CARP_DIR("GlVTKConvert"); };
        inline QString GlRuleFibres() { return CARP_DIR("GlRuleFibres"); };
        inline QString GlElemCenters() { return CARP_DIR("GlElemCenters"); };
        inline QString carp_pt() { return CARP_DIR("carp.pt"); };
        inline QString igbextract() { return CARP_DIR("igbextract"); };

    protected:
        // Segmentation Utilities
        void GetLabels(mitk::Image::Pointer seg, std::vector<int>& labels, int background=0);
        mitk::Image::Pointer ExtractSingleLabel(mitk::Image::Pointer seg, int label, bool binarise=true);
        mitk::Image::Pointer BwLabelN(mitk::Image::Pointer seg, std::vector<int>& labels);

    private:
        QString _carp_dir = "";
};
#endif // CemrgFourChamberTools_h