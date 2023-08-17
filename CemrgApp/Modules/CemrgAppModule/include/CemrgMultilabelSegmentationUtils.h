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

#ifndef CemrgMultilabelSegmentationUtils_h
#define CemrgMultilabelSegmentationUtils_h

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

class MITKCEMRGAPPMODULE_EXPORT CemrgMultilabelSegmentationUtils {
    typedef itk::Image<uint8_t, 3> ImageType;
    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    
    public:
        CemrgMultilabelSegmentationUtils();
        ~CemrgMultilabelSegmentationUtils();

        inline void SetDebugDir(QString dir) { debugDir = dir; };
        inline QString GetDebugDir() { return debugDir; };

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

    protected:

    private : 
        QString debugDir;
       
};
#endif // CemrgMultilabelSegmentationUtils_h
