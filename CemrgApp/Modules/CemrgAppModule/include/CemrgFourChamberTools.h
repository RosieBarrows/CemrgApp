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

        inline void SetDebug(bool debug, std::string debugDirectory) { debug = debug; debugDir = debugDirectory; };
        inline void SetDebugOn(std::string debugDirectory) { SetDebug(true, debugDirectory); };

        inline void SetDebugDir(std::string debugDirectory) { debugDir = debugDirectory; };
        inline std::string GetDebugDir() { return debugDir; };
        inline QString QGetDebugDir() { return QString::fromStdString(debugDir); };

        inline void UpdateStep(mitk::Image::Pointer image, SegmentationStep step) { currentImage = image; currentStep = step; };
        inline void SetLabel(LabelsType stt, int label) { chosenLabels.Set(stt, label); };

        // Segmentation Utilities
        bool CheckExisting(mitk::Image::Pointer seg, int queryLabel);
        void ExploreLabelsToSplit(mitk::Image::Pointer seg, std::vector<int>& labels); 
        mitk::Image::Pointer SplitLabelsOnRepeat(mitk::Image::Pointer seg, int label, unsigned int radius=3);
        
        mitk::Image::Pointer ReplaceLabel(mitk::Image::Pointer seg, int oldLabel, int newLabel);
        inline mitk::Image::Pointer RemoveLabel(mitk::Image::Pointer seg, int label) { return ReplaceLabel(seg, label, 0); };

        enum MaskLabelBehaviour { ZEROS, ONLY, REPLACE, FORBID };
        mitk::Image::Pointer AddMaskToSegmentation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel, MaskLabelBehaviour mlb, std::vector<int> labelsToProcess = std::vector<int>());
        inline mitk::Image::Pointer AddMask(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel) { return AddMaskToSegmentation(seg, mask, newLabel, MaskLabelBehaviour::ZEROS); };
        inline mitk::Image::Pointer AddMaskReplace(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel) { return AddMaskToSegmentation(seg, mask, newLabel, MaskLabelBehaviour::REPLACE); };
        inline mitk::Image::Pointer AddMaskReplaceOnly(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel, std::vector<int> labelsToProcess) { return AddMaskToSegmentation(seg, mask, newLabel, MaskLabelBehaviour::REPLACE, labelsToProcess); };
        inline mitk::Image::Pointer AddMaskReplaceExcept(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel, std::vector<int> labelsToProcess) { return AddMaskToSegmentation(seg, mask, newLabel, MaskLabelBehaviour::FORBID, labelsToProcess); };

        mitk::Image::Pointer LabelMaskAndOperation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int oldLabel, int newLabel);
        mitk::Image::Pointer ConnectedComponent(mitk::Image::Pointer seg, std::vector<unsigned int> seedIdx, int label, bool keep=false);
        inline mitk::Image::Pointer ConnectedComponentKeep(mitk::Image::Pointer seg, std::vector<unsigned int> seedIdx, int label) { return ConnectedComponent(seg, seedIdx, label, true); };

        // helper functions
        void GetLabels(mitk::Image::Pointer seg, std::vector<int>& labels, int background=0);
        mitk::Image::Pointer ExtractSingleLabel(mitk::Image::Pointer seg, int label, bool binarise=true);
        mitk::Image::Pointer BwLabelN(mitk::Image::Pointer seg, std::vector<int>& labels);

        bool GetLabelCentreOfMassIndex(mitk::Image::Pointer seg, int label, std::vector<unsigned int> &cogIndx);
        bool GetLabelCentreOfMass(mitk::Image::Pointer seg, int label, std::vector<double> &cog);

        void WorldToIndex(mitk::Image::Pointer image, std::vector<double> world, std::vector<unsigned int>& index);
        void IndexToWorld(mitk::Image::Pointer image, std::vector<unsigned int> index, std::vector<double> &world);

        // Manipulation utilities
        void SetCylinders(QJsonObject json);
        void SetSlicers(QJsonObject json);
        void SetValvePoints(QJsonObject json);

        inline bool CylindersSet() { return _cylinders.IsPointSet(); };
        inline bool SlicersSet() { return _slicers.IsPointSet(); };
        inline bool ValvePointsSet() { return _valvePoints.IsPointSet(); };

        mitk::Image::Pointer Cylinder(mitk::Image::Pointer seg, QString ptPrefix, double slicerRadius, double slicerHeight, ManualPoints mpl, QString saveAs = "");
        mitk::Image::Pointer CreateSvcIvc(std::vector<mitk::Image::Pointer> images, int RspvLabel=10, int SvcLabel=13, int IvcLabel=14);
        mitk::Image::Pointer CropSvcIvc(std::vector<mitk::Image::Pointer> images,
                                        std::vector<unsigned int> seedSVC,
                                        std::vector<unsigned int> seedIVC,
                                        int aortaSlicerLabel,
                                        int PArtSlicerLabel);

    protected:
        mitk::Image::Pointer S2B(std::vector<unsigned int> seedSVC);
        mitk::Image::Pointer S2C(std::vector<unsigned int> seedIVC);
        mitk::Image::Pointer S2D(mitk::Image::Pointer aorta, mitk::Image::Pointer PArt, mitk::Image::Pointer svc_slicer, mitk::Image::Pointer ivc_slicer, int aortaSlicerLabel, int PArtSlicerLabel);
        mitk::Image::Pointer S2E(std::vector<unsigned int> seedSVC);
        mitk::Image::Pointer S2F(std::vector<unsigned int> seedIVC);

    private : 
        CylinderPointsType _cylinders;
        SlicersPointsType _slicers;
        ValvePlainsPointsType _valvePoints;

        bool debug;
        std::string debugDir;

        mitk::Image::Pointer currentImage;
        SegmentationStep currentStep;
        LabelsStruct chosenLabels;
};
#endif // CemrgFourChamberTools_h
