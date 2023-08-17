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
#include <itkConnectedThresholdImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>

#include "CemrgCommonUtils.h"
#include "CemrgAtrialTools.h"
#include "CemrgMultilabelSegmentationUtils.h"

CemrgMultilabelSegmentationUtils::CemrgMultilabelSegmentationUtils(){
    
}

CemrgMultilabelSegmentationUtils::~CemrgMultilabelSegmentationUtils(){

}


// Segmentation Utilities
void CemrgMultilabelSegmentationUtils::ExploreLabelsToSplit(mitk::Image::Pointer seg, std::vector<int> &labels) {
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

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::SplitLabelsOnRepeat(mitk::Image::Pointer seg, int label, unsigned int radius) {
    
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

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::ReplaceLabel(mitk::Image::Pointer seg, int oldLabel, int newLabel) {
    ImageType::Pointer itkImage = ImageType::New();
    CastToItkImage(seg, itkImage);

    IteratorType segIt(itkImage, itkImage->GetLargestPossibleRegion());
    for (segIt.GoToBegin(); !segIt.IsAtEnd(); ++segIt) {
        if (segIt.Get() == oldLabel) {
            segIt.Set(newLabel);
        }
    }
    
    mitk::Image::Pointer newSeg = mitk::ImportItkImage(itkImage)->Clone();
    newSeg->SetGeometry(seg->GetGeometry());
    
    return newSeg;
}

bool CemrgMultilabelSegmentationUtils::CheckExisting(mitk::Image::Pointer seg, int queryLabel) {
    std::vector<int> labelsInSeg;
    GetLabels(seg, labelsInSeg);

    MITK_INFO << ("Checking if label " + QString::number(queryLabel) + " exists in segmentation").toStdString();
    auto it = std::find(labelsInSeg.begin(), labelsInSeg.end(), queryLabel);

    return (it != labelsInSeg.end()); 
}

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::AddMaskToSegmentation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel, MaskLabelBehaviour mlb, std::vector<int> labelsToProcess) {

    ImageType::Pointer itkImage = ImageType::New();
    ImageType::Pointer itkMask = ImageType::New();
    CastToItkImage(seg, itkImage);
    CastToItkImage(mask, itkMask);

    IteratorType segIt(itkImage, itkImage->GetLargestPossibleRegion());
    IteratorType maskIt(itkMask, itkMask->GetLargestPossibleRegion());

    for (segIt.GoToBegin(), maskIt.GoToBegin(); !segIt.IsAtEnd(); ++segIt, ++maskIt) {
        bool conditionToProcess = false;
        uint8_t value = segIt.Get();

        if (maskIt.Get() > 0) {
            auto testValueInLabels = std::find(labelsToProcess.begin(), labelsToProcess.end(), value);
            switch(mlb) { 
                case MaskLabelBehaviour::ZEROS: 
                    conditionToProcess = (value == 0);
                    break;
                case MaskLabelBehaviour::ONLY:
                    conditionToProcess = (value == 0) || (testValueInLabels != labelsToProcess.end());
                    break;
                case MaskLabelBehaviour::REPLACE:
                    conditionToProcess = true;
                    break;
                case MaskLabelBehaviour::FORBID:
                    conditionToProcess = (testValueInLabels == labelsToProcess.end());
                    break;
                default:
                    conditionToProcess = false;
                }

            if (conditionToProcess) {
                segIt.Set(newLabel);
            }
        }
    }

    mitk::Image::Pointer newSeg = mitk::ImportItkImage(itkImage)->Clone();
    newSeg->SetGeometry(seg->GetGeometry());
    return newSeg;

}

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::LabelMaskAndOperation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int oldLabel, int newLabel) {
    ImageType::Pointer itkImage = ImageType::New();
    ImageType::Pointer itkMask = ImageType::New();
    CastToItkImage(seg, itkImage);
    CastToItkImage(mask, itkMask);

    IteratorType segIt(itkImage, itkImage->GetLargestPossibleRegion());
    IteratorType maskIt(itkMask, itkMask->GetLargestPossibleRegion());

    for (segIt.GoToBegin(), maskIt.GoToBegin(); !segIt.IsAtEnd(); ++segIt, ++maskIt) {
        if (maskIt.Get() > 0) {
            int valueToSet = (segIt.Get() == oldLabel) ? newLabel : 0;
            segIt.Set(valueToSet);
        }
    }

    mitk::Image::Pointer newSeg = mitk::ImportItkImage(itkImage)->Clone();
    newSeg->SetGeometry(seg->GetGeometry());
    return newSeg;
}

void CemrgMultilabelSegmentationUtils::GetLabels(mitk::Image::Pointer seg, std::vector<int> &labels, int background) {
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

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::ExtractSingleLabel(mitk::Image::Pointer seg, int label, bool binarise) {
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

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::BwLabelN(mitk::Image::Pointer seg, std::vector<int> &labels) {
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

bool CemrgMultilabelSegmentationUtils::GetLabelCentreOfMassIndex(mitk::Image::Pointer seg, int label, std::vector<unsigned int> &cogIndx) {
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
        cogIndx.push_back((bb[2*ix] + bb[2*ix + 1]) / 2);
    }

    return true;
}

bool CemrgMultilabelSegmentationUtils::GetLabelCentreOfMass(mitk::Image::Pointer seg, int label, std::vector<double> &cog) {
    std::vector<unsigned int> cogIndx;
    bool result = GetLabelCentreOfMassIndex(seg, label, cogIndx);    
    if (!result) {
        return false;
    }

    IndexToWorld(seg, cogIndx, cog);

    return true;
}

void CemrgMultilabelSegmentationUtils::WorldToIndex(mitk::Image::Pointer image, std::vector<double> world, std::vector<unsigned int>& index) {
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    ImageType::PointType point;
    ImageType::IndexType indexType;
    for (unsigned int ix = 0; ix < 3; ix++) {
        point[ix] = world.at(ix);
    }

    itkImage->TransformPhysicalPointToIndex(point, indexType);

    for (unsigned int ix = 0; ix < 3; ix++) {
        index.push_back(indexType[ix]);
    }
}

void CemrgMultilabelSegmentationUtils::IndexToWorld(mitk::Image::Pointer image, std::vector<unsigned int> index, std::vector<double> &world) {
    
    for (unsigned int ix = 0; ix < 3; ix++) {
        double spacing = image->GetGeometry()->GetSpacing()[ix];
        double origin = image->GetGeometry()->GetOrigin()[ix];

        world.push_back( (index[ix] * spacing) + origin );   
    }

}

mitk::Image::Pointer CemrgMultilabelSegmentationUtils::ConnectedComponent(mitk::Image::Pointer seg, std::vector<unsigned int> seedIdx, int layer, bool keep) {
    using ConnectedThresholdType = itk::ConnectedThresholdImageFilter<ImageType, ImageType>;
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(seg, itkImage);
    ImageType::IndexType seedPointIndex;
    seedPointIndex[0] = seedIdx.at(0);
    seedPointIndex[1] = seedIdx.at(1);
    seedPointIndex[2] = seedIdx.at(2);
    
    ConnectedThresholdType::Pointer cc = ConnectedThresholdType::New();
    cc->SetInput(itkImage);
    cc->SetSeed(seedPointIndex);
    cc->SetLower(layer);
    cc->SetUpper(layer);

    mitk::Image::Pointer ccImage = mitk::ImportItkImage(cc->GetOutput())->Clone();
    ccImage->SetGeometry(seg->GetGeometry());

    
    MITK_INFO << "Saving connected threshold.";
    mitk::IOUtil::Save(ccImage, (debugDir + "/CC.nii").toStdString());

    int replaceValue = 0;
    if (keep) {
        seg = RemoveLabel(seg, layer);
        replaceValue = layer;
    }

    return AddMaskReplace(seg, ccImage, replaceValue);

}
