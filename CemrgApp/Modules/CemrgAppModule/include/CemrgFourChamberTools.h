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
#include <QJsonObject>

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

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberTools {
    typedef itk::Image<uint8_t, 3> ImageType;
    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    
    public:
        CemrgFourChamberTools();
        ~CemrgFourChamberTools();

        // Segmentation Utilities
        void ExploreLabelsToSplit(mitk::Image::Pointer seg, std::vector<int>& labels); 
        mitk::Image::Pointer SplitLabelsOnRepeat(mitk::Image::Pointer seg, int label, unsigned int radius=3);
        mitk::Image::Pointer SwapLabel(mitk::Image::Pointer seg, int label1, int label2, bool checkExisting = true, bool overwrite = false);
        mitk::Image::Pointer AddMaskToSegmentation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int label, std::vector<int> labelsToIgnore = std::vector<int>());
        inline mitk::Image::Pointer RemoveLabel(mitk::Image::Pointer seg, int label) { return SwapLabel(seg, label, 0, false, true); };

        void GetLabels(mitk::Image::Pointer seg, std::vector<int>& labels, int background=0);
        mitk::Image::Pointer ExtractSingleLabel(mitk::Image::Pointer seg, int label, bool binarise=true);
        mitk::Image::Pointer BwLabelN(mitk::Image::Pointer seg, std::vector<int>& labels);

        bool GetLabelCentreOfMassIndex(mitk::Image::Pointer seg, int label, std::vector<unsigned int> &cogIndx);
        bool GetLabelCentreOfMass(mitk::Image::Pointer seg, int label, std::vector<double> &cog);

        bool WorldToIndex(mitk::Image::Pointer image, std::vector<double> world, std::vector<unsigned int>& index);
        bool IndexToWorld(mitk::Image::Pointer image, std::vector<unsigned int> index, std::vector<double> &world);

        // Manipulation utilities
        void SetCylinders(QJsonObject json);
        void SetSlicers(QJsonObject json);
        void SetValvePoints(QJsonObject json);

        inline bool CylindersSet() { return _cylinders.pointsSet; };
        inline bool SlicersSet() { return _slicers.pointsSet; };
        inline bool ValvePointsSet() { return _valvePoints.pointsSet; };

        mitk::Image::Pointer Cylinder(mitk::Image::Pointer seg, QString ptPrefix, double slicerRadius, double slicerHeight, QString saveAs = "");
        mitk::Image::Pointer CreateSvcIvc(std::vector<mitk::Image::Pointer> images, int RspvLabel=10, int SvcLabel=13, int IvcLabel=14);
    protected:
        

    private:
        CylinderPointsType _cylinders;
        SlicersPointsType _slicers;
        ValvePlainsPointsType _valvePoints;
};
#endif // CemrgFourChamberTools_h
