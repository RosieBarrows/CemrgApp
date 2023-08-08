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
#include <mitkImagePixelWriteAccessor.h>

// ITK
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>

#include "CemrgCommonUtils.h"
#include "CemrgAtrialTools.h"
#include "CemrgFourChamberTools.h"

CemrgFourChamberTools::CemrgFourChamberTools(){

}

CemrgFourChamberTools::~CemrgFourChamberTools(){

}


// Segmentation Utilities
void CemrgFourChamberTools::ExploreLabelsToSplit(mitk::Image::Pointer seg, std::vector<int> &labels) {
    std::vector<int> labelsInSeg;
    GetLabels(seg, labelsInSeg);

    for (long unsigned int ix = 0; ix < labelsInSeg.size(); ix++) { 
        int label = labelsInSeg.at(ix);
        std::vector<int> tagsInLabel;
        mitk::Image::Pointer bwlabelnInLabel = BwLabelN(ExtractSingleLabel(seg, label), tagsInLabel);
        if (tagsInLabel.size() > 1) {
            labels.push_back(label);
        }
    }
}

mitk::Image::Pointer CemrgFourChamberTools::SplitLabelsOnRepeat(mitk::Image::Pointer seg, int label, unsigned int radius) {
    
    std::vector<int> labelsInSeg;
    GetLabels(seg, labelsInSeg);
    // remove label from labelsInSeg
    labelsInSeg.erase(std::remove(labelsInSeg.begin(), labelsInSeg.end(), label), labelsInSeg.end());

    std::unique_ptr<CemrgAtrialTools> imatools = std::unique_ptr<CemrgAtrialTools>(new CemrgAtrialTools());
    imatools->SetDebugModeOff();

    mitk::Image::Pointer imageLabel = ExtractSingleLabel(seg, label, true);
    ImageType::Pointer itkLabel = ImageType::New();
    CastToItkImage(imageLabel, itkLabel);

    using StrElType = itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>;
    using OpeningType = itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType, StrElType>;
    StrElType strel;
    strel.SetRadius(radius);
    strel.CreateStructuringElement();
    OpeningType::Pointer imopen = OpeningType::New();
    imopen->SetInput(itkLabel);
    imopen->SetKernel(strel);
    imopen->Update();

    std::vector<int> tagsInLabel;
    mitk::Image::Pointer ccImLabel = BwLabelN(mitk::ImportItkImage(imopen->GetOutput()) , tagsInLabel);

    if (tagsInLabel.size() == 1) {
        return seg;
    }

    ImageType::Pointer itkImage = ImageType::New();
    CastToItkImage(seg, itkImage);
    ImageType::Pointer itkCcImLabel = ImageType::New();
    CastToItkImage(ccImLabel, itkCcImLabel);

    // Remove label from seg
    IteratorType it(itkImage, itkImage->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        if (it.Get() == label) {
            it.Set(0);
        }
    }

    for (long unsigned int ix = 0; ix < tagsInLabel.size(); ix++) {
        int tag = tagsInLabel.at(ix);
        int newLabel = (ix==0) ? label : label*10 + (tag-1);

        int qx = 1;
        // Use std::find to check if the variable is in the vector
        auto it = std::find(tagsInLabel.begin(), tagsInLabel.end(), newLabel);
        while (it != tagsInLabel.end()) {
            newLabel = label*std::pow(10, qx) + (tag-1);
            qx++;
            it = std::find(tagsInLabel.begin(), tagsInLabel.end(), newLabel);
        }

        IteratorType ccIt(itkCcImLabel, itkCcImLabel->GetLargestPossibleRegion());
        IteratorType segIt(itkImage, itkImage->GetLargestPossibleRegion());

        segIt.GoToBegin();
        for (ccIt.GoToBegin(); !ccIt.IsAtEnd(); ++ccIt) {
            if (ccIt.Get() == tag) {
                segIt.Set(newLabel);
            }
            ++segIt;
        }
    }

    mitk::Image::Pointer newSeg = mitk::ImportItkImage(itkImage)->Clone();
    return newSeg;

}

mitk::Image::Pointer CemrgFourChamberTools::SwapLabel(mitk::Image::Pointer seg, int label1, int label2, bool checkExisting, bool overwrite) {
    std::vector<int> labelsInSeg;
    GetLabels(seg, labelsInSeg);

    if (checkExisting) {
        MITK_INFO << ("Checking if label " + QString::number(label2) + " exists in segmentation").toStdString();
        auto it = std::find(labelsInSeg.begin(), labelsInSeg.end(), QString::number(label2));
        if (it != labelsInSeg.end() && !overwrite) {
            MITK_INFO << ("Label " + QString::number(label2) + " already exists in segmentation").toStdString();
            return seg;
        }
    }

    MITK_INFO(overwrite) << ("Overwritting label " + QString::number(label2)).toStdString();
    ImageType::Pointer itkImage = ImageType::New();
    CastToItkImage(seg, itkImage);

    IteratorType segIt(itkImage, itkImage->GetLargestPossibleRegion());
    for (segIt.GoToBegin(); !segIt.IsAtEnd(); ++segIt) {
        if (segIt.Get() == label1) {
            segIt.Set(label2);
        } else if (segIt.Get() == label2) {
            segIt.Set(label1);
        }
    }
    
    mitk::Image::Pointer newSeg = mitk::ImportItkImage(itkImage)->Clone();
    newSeg->SetGeometry(seg->GetGeometry());
    return newSeg;
    
}

void CemrgFourChamberTools::GetLabels(mitk::Image::Pointer seg, std::vector<int> &labels, int background) {
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(seg, itkImage);

    IteratorType it(itkImage, itkImage->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        int label = it.Get();
        if (label != background) {
            labels.push_back(label);
        }
    }

    std::sort(labels.begin(), labels.end());
    labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
}

mitk::Image::Pointer CemrgFourChamberTools::ExtractSingleLabel(mitk::Image::Pointer seg, int label, bool binarise) {
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(seg, itkImage);

    ImageType::Pointer outImage = ImageType::New();
    outImage->SetRegions(itkImage->GetLargestPossibleRegion());
    outImage->SetSpacing(itkImage->GetSpacing());
    outImage->SetOrigin(itkImage->GetOrigin());
    outImage->Allocate();
    outImage->FillBuffer(0);

    IteratorType it(itkImage, itkImage->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        int valueAtIter = it.Get();
        int valueInOutput = binarise ? 1 : valueAtIter;
        if (valueAtIter == label) {
            outImage->SetPixel(it.GetIndex(), valueInOutput);
        }
    }

    mitk::Image::Pointer outSeg = mitk::ImportItkImage(outImage)->Clone();

    return outSeg;
}

mitk::Image::Pointer CemrgFourChamberTools::BwLabelN(mitk::Image::Pointer seg, std::vector<int> &labels) {
    // assumes seg is a binary image
    using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<ImageType, ImageType>;
    using RelabelFilterType = itk::RelabelComponentImageFilter<ImageType, ImageType> ;

    mitk::MorphologicalOperations::Opening(seg, 1, mitk::MorphologicalOperations::StructuralElementType::Ball);

    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(seg, itkImage);

    ConnectedComponentImageFilterType::Pointer conn1 = ConnectedComponentImageFilterType::New();
    conn1->SetInput(itkImage);
    conn1->Update();

    RelabelFilterType::Pointer labelled = RelabelFilterType::New();
    labelled->SetInput(conn1->GetOutput());
    labelled->Update();

    mitk::Image::Pointer outImage = mitk::ImportItkImage(labelled->GetOutput())->Clone();
    
    GetLabels(outImage, labels);
    
    return outImage;
}

bool CemrgFourChamberTools::GetLabelCentreOfMassIndex(mitk::Image::Pointer seg, int label, std::vector<unsigned int> &cogIndx) {
    using LabelStatisticsFilterType = itk::LabelStatisticsImageFilter<ImageType, ImageType>;

    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(seg, itkImage);

    LabelStatisticsFilterType::Pointer labelStats = LabelStatisticsFilterType::New();
    labelStats->SetInput(itkImage);
    labelStats->SetLabelInput(itkImage); // Use the segmentation image as label
    labelStats->Update();

    LabelStatisticsFilterType::LabelPixelType labelValue = label;
    LabelStatisticsFilterType::BoundingBoxType bb = labelStats->GetBoundingBox(labelValue);

    if (bb.empty()) {
        return false;
    }
    
    for (unsigned int ix = 0; ix < 3; ix++) {
        double spacing = itkImage->GetSpacing()[ix];
        double origin = itkImage->GetOrigin()[ix];
        cogIndx.push_back((bb[2*ix] + bb[2*ix + 1]) / 2);
    }

    return true;
}

bool CemrgFourChamberTools::GetLabelCentreOfMass(mitk::Image::Pointer seg, int label, std::vector<double> &cog) {
    std::vector<unsigned int> cogIndx;
    bool result = GetLabelCentreOfMassIndex(seg, label, cogIndx);    
    if (!result) {
        return false;
    }

    for (unsigned int ix = 0; ix < 3; ix++) {
        double spacing = seg->GetGeometry()->GetSpacing()[ix];
        double origin = seg->GetGeometry()->GetOrigin()[ix];

        cog.push_back( (cogIndx[ix] * spacing) + origin );   
    }

    std::cout << "Centre of mass: " << cog[0] << ", " << cog[1] << ", " << cog[2] << '\n';

    return true;
}

bool CemrgFourChamberTools::WorldToIndex(mitk::Image::Pointer image, std::vector<double> world, std::vector<unsigned int>& index) {
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    return false;
}

bool CemrgFourChamberTools::IndexToWorld(mitk::Image::Pointer image, std::vector<unsigned int> index, std::vector<double> &world) {
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    return false;
}
