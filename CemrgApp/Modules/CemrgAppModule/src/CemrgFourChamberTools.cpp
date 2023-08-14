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
#include "CemrgFourChamberTools.h"

CemrgFourChamberTools::CemrgFourChamberTools(){
    _cylinders = CylinderPointsType();
    _slicers = SlicersPointsType();
    _valvePoints = ValvePlainsPointsType();
    debug = false;
    debugDir = "";
    currentStep = SegmentationStep::NONE;
    currentImage = nullptr;
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

mitk::Image::Pointer CemrgFourChamberTools::ReplaceLabel(mitk::Image::Pointer seg, int oldLabel, int newLabel) {
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



bool CemrgFourChamberTools::CheckExisting(mitk::Image::Pointer seg, int queryLabel) {
    std::vector<int> labelsInSeg;
    GetLabels(seg, labelsInSeg);

    MITK_INFO << ("Checking if label " + QString::number(queryLabel) + " exists in segmentation").toStdString();
    auto it = std::find(labelsInSeg.begin(), labelsInSeg.end(), queryLabel);

    return (it != labelsInSeg.end()); 
}

mitk::Image::Pointer CemrgFourChamberTools::AddMaskToSegmentation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int newLabel, MaskLabelBehaviour mlb, std::vector<int> labelsToProcess) {

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

mitk::Image::Pointer CemrgFourChamberTools::LabelMaskAndOperation(mitk::Image::Pointer seg, mitk::Image::Pointer mask, int oldLabel, int newLabel) {
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

    IndexToWorld(seg, cogIndx, cog);

    return true;
}

void CemrgFourChamberTools::WorldToIndex(mitk::Image::Pointer image, std::vector<double> world, std::vector<unsigned int>& index) {
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

void CemrgFourChamberTools::IndexToWorld(mitk::Image::Pointer image, std::vector<unsigned int> index, std::vector<double> &world) {
    
    for (unsigned int ix = 0; ix < 3; ix++) {
        double spacing = image->GetGeometry()->GetSpacing()[ix];
        double origin = image->GetGeometry()->GetOrigin()[ix];

        world.push_back( (index[ix] * spacing) + origin );   
    }

}

void CemrgFourChamberTools::SetCylinders(QJsonObject json) {
    QStringList keys = json.keys();
    foreach (QString key, keys) {
        _cylinders.SetPoint(json, key);
    }
    _cylinders.PointSetOn();
}

void CemrgFourChamberTools::SetSlicers(QJsonObject json) {
    QStringList keys = json.keys();
    foreach (QString key, keys) {
        _slicers.SetPoint(json, key);
    }
    _slicers.PointSetOn();
}

void CemrgFourChamberTools::SetValvePoints(QJsonObject json) {
    QStringList keys = json.keys();
    foreach (QString key, keys) {
        _valvePoints.SetPoint(json, key);
    }
    _valvePoints.PointSetOn();
}

struct PointList {
    std::vector<double> pt1, pt2, pt3, v1, v2, _normal;
    PointList() 
        : pt1(3, 0.0), 
          pt2(3, 0.0), 
          pt3(3, 0.0), 
          v1(3, 0.0), v2(3, 0.0), _normal(3, 0.0) {};
    
    void FillPoint(int whichPoint, std::vector<double> point) {
        switch (whichPoint) {
            case 1: pt1 = point; break;
            case 2: pt2 = point; break;
            case 3: pt3 = point; break;
            default: break;
        }
    }

    std::vector<double> pt(int whichPoint) {
        switch (whichPoint) {
            case 1: return pt1;
            case 2: return pt2;
            case 3: return pt3;
            default: return std::vector<double>(3, 0.0);
        }
    }

    void normalise(std::vector<double>& v) {
        double norm = std::sqrt(std::pow(v.at(0), 2) + std::pow(v.at(1), 2) + std::pow(v.at(2), 2));
        for (int ix = 0; ix < 3; ix++) {
            v.at(ix) /= norm;
        }
    }

    std::vector<double> v(int a, int b) {
        std::vector<double> v(3, 0.0);
        for (int ix = 0; ix < 3; ix++) {
            v.at(ix) = pt(b).at(ix) - pt(a).at(ix);
        }
        normalise(v);

        return v;
    }

    void cross() {
        _normal.at(0) = v1.at(1)*v2.at(2) - v1.at(2)*v2.at(1);
        _normal.at(1) = v1.at(2)*v2.at(0) - v1.at(0)*v2.at(2);
        _normal.at(2) = v1.at(0)*v2.at(1) - v1.at(1)*v2.at(0);
        normalise(_normal);
    }

    void SetVectors() {
        v1 = v(1, 2);
        v2 = v(1, 3);
    }

    std::vector<double> normal() { return _normal; };

    std::vector<double> getCross(std::vector<double> a, std::vector<double> b) {
        std::vector<double> cross(3, 0.0);
        cross.at(0) = a.at(1)*b.at(2) - a.at(2)*b.at(1);
        cross.at(1) = a.at(2)*b.at(0) - a.at(0)*b.at(2);
        cross.at(2) = a.at(0)*b.at(1) - a.at(1)*b.at(0);
        return cross;
    }

    double getNorm(std::vector<double> a) {
        return std::sqrt(std::pow(a.at(0), 2) + std::pow(a.at(1), 2) + std::pow(a.at(2), 2));
    }
};

/// @brief 
/// @param seg - input segmentation
/// @param ptPrefix - prefix of the points to be used CYLINDERS(SVC, IVC, ...) SLICERS(SVC_slicer, IVC_Slicer)
/// @param slicerRadius 
/// @param slicerHeight 
/// @param mpl - ManualPoints
/// @param saveAs 
/// @return 
mitk::Image::Pointer CemrgFourChamberTools::Cylinder(mitk::Image::Pointer seg, QString ptPrefix, double slicerRadius, double slicerHeight, ManualPoints mpl, QString saveAs) {

    // Get keys whether it is cylinders, slicers, or valve_plains
    ManualPointsStruct segIds;
    QStringList keys = segIds.GetPointLabelOptions(mpl);

    BasePointsType* bpt;
    switch (mpl) { 
        case ManualPoints::CYLINDERS: 
            bpt = &_cylinders; 
            break;
        case ManualPoints::SLICERS:
            bpt = &_slicers;
            break;
        case ManualPoints::VALVE_PLAINS:
            bpt = &_valvePoints;
            break;
        default: 
            MITK_ERROR << "Invalid ManualPoints";
            return nullptr;
    }

    if (!keys.contains(ptPrefix)) {
        MITK_WARN << ("Prefix " + ptPrefix + " not found in segmentation").toStdString();
        return nullptr;
    }
    
    // remove keys that do not contain the prefix
    for (int ix = 0; ix < keys.size(); ix++) {
        if (!keys.at(ix).contains(ptPrefix)) {
            keys.removeAt(ix);
            ix--;
        }
    }

    if (keys.size() != 3) {
        MITK_WARN << ("Prefix " + ptPrefix + " does not have 3 points").toStdString();
        return nullptr;
    }

    double origin[3], spacing[3];
    seg->GetGeometry()->GetOrigin().ToArray(origin);
    seg->GetGeometry()->GetSpacing().ToArray(spacing);

    ImageType::Pointer outputImg = ImageType::New();
    mitk::CastToItkImage(seg->Clone(), outputImg);

    using IteratorType = itk::ImageRegionIterator<ImageType>;
    IteratorType imIter(outputImg, outputImg->GetLargestPossibleRegion());
    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        imIter.Set(0);
        ++imIter;
    }

    int numPoints = keys.size(); 

    std::vector<double> cog(3, 0.0);
    PointList pts;
    int pointId = 1;

    foreach (QString key, keys) {

        PointsNamesType pointType = bpt->FromKey(key);

        std::vector<double> thisPoint(3);
        for (int ix = 0; ix < 3; ix++) {
            thisPoint.at(ix) = bpt->GetPointAt(pointType, ix);
            cog.at(ix) += thisPoint.at(ix);
        }
        pts.FillPoint(pointId, thisPoint);
        pointId++;
    }

    for (int ix = 0; ix < 3; ix++) {
        cog.at(ix) /= numPoints;
    }

    pts.SetVectors();
    pts.cross();

    std::vector<double> p1(3), p2(3), normalVector(3);
    double n_x, n_y, n_z;

    for (int ix = 0; ix < 3; ix++) {
        p1.at(ix) = cog.at(ix) - pts.normal().at(ix) * (slicerHeight/2); 
        p2.at(ix) = cog.at(ix) + pts.normal().at(ix) * (slicerHeight/2); 

        normalVector.at(ix) = p2.at(ix) - p1.at(ix);
    }

    n_x = seg->GetDimensions()[0];
    n_y = seg->GetDimensions()[1];
    n_z = seg->GetDimensions()[2];

    double cubeSize = std::max(slicerHeight, slicerRadius);
    cubeSize += (slicerHeight > slicerRadius) ? 10 : 30;

    std::vector<int> xCubeCoord, yCubeCoord, zCubeCoord;

    for (int dimension = 0; dimension < 3;dimension++) {
        int upperBound = (dimension == 0) ? n_x : ((dimension == 1) ? n_y : n_z);

        for (int ix = 0; ix < upperBound; ix++) {
            double coord = origin[dimension] + ix * spacing[dimension];
            if (std::abs(coord - cog.at(dimension)) < (cubeSize / 2)) {
                switch (dimension) {
                    case 0: xCubeCoord.push_back(ix); break;
                    case 1: yCubeCoord.push_back(ix); break;
                    case 2: zCubeCoord.push_back(ix); break;
                    default: break;
                }
            }
        }
    }

    foreach (int xcoord, xCubeCoord) {
        foreach (int ycoord, yCubeCoord) { 
            foreach (int zcoord, zCubeCoord) {
                std::vector<double> testPts(3), v1(3), v2(3);
                for (int ix = 0; ix<3; ix++) {
                    testPts.at(ix) = origin[ix] + xcoord * spacing[ix];
                    v1.at(ix) = testPts.at(ix) - p1.at(ix);
                    v2.at(ix) = testPts.at(ix) - p2.at(ix);
                }
                double dot1, dot2;
                dot1 = v1.at(0)*normalVector.at(0) + v1.at(1)*normalVector.at(1) + v1.at(2)*normalVector.at(2);
                dot2 = v2.at(0)*normalVector.at(0) + v2.at(1)*normalVector.at(1) + v2.at(2)*normalVector.at(2);
                if (dot1 >= 0 && dot2 <= 0) {
                    pts.normalise(normalVector);
                    double testRadius = pts.getNorm(pts.getCross(testPts, normalVector));
                    if (testRadius <= slicerRadius) {
                        ImageType::IndexType index;
                        index[0] = xcoord;
                        index[1] = ycoord;
                        index[2] = zcoord;
                        outputImg->SetPixel(index, 1);
                    }
                }
            }
        }
    }

    mitk::Image::Pointer output = mitk::ImportItkImage(outputImg)->Clone();
    output->SetGeometry(seg->GetGeometry());

    if (output) {
        if (!saveAs.isEmpty()) {
            mitk::IOUtil::Save(output, saveAs.toStdString());
        }
    } else {
        MITK_WARN << "Could not create cylinder: ";
        MITK_WARN(!saveAs.isEmpty()) << saveAs.toStdString();
    }

    return output;
}

mitk::Image::Pointer CemrgFourChamberTools::CreateSvcIvc(std::vector<mitk::Image::Pointer> images, int RspvLabel, int SvcLabel, int IvcLabel) {
    if (images.size() != 3) {
        MITK_WARN << "Need 3 images to create SVC/IVC";
        return nullptr;
    }

    mitk::Image::Pointer seg = images.at(0);
    mitk::Image::Pointer svc = images.at(1);
    mitk::Image::Pointer ivc = images.at(2);

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(RspvLabel);
    mitk::Image::Pointer svcIvc = AddMaskReplaceOnly(seg, svc, SvcLabel, labelsToProcess);

    if (svcIvc == nullptr) {
        return nullptr;
    }

    return AddMask(svcIvc, ivc, IvcLabel);
}

mitk::Image::Pointer CemrgFourChamberTools::S2B(std::vector<unsigned int> seedSVC) {
    if (currentStep != SegmentationStep::S2A) {
        MITK_WARN << "Need to run S2A first";
        return nullptr;
    }

    mitk::Image::Pointer s2b = ConnectedComponentKeep(currentImage, seedSVC, chosenLabels.Get(LabelsType::SVC));
    if (debug) {
        MITK_INFO << "Saving connected threshold.";
        mitk::IOUtil::Save(s2b, debugDir + "/seg_s2b.nii");
    }

    UpdateStep(s2b, SegmentationStep::S2B);

    return s2b;
}

mitk::Image::Pointer CemrgFourChamberTools::S2C(std::vector<unsigned int> seedIVC) {
    if (currentStep != SegmentationStep::S2B) {
        MITK_WARN << "Need to run S2B first";
        return nullptr;
    }

    mitk::Image::Pointer s2c = ConnectedComponentKeep(currentImage, seedIVC, chosenLabels.Get(LabelsType::IVC));
    if (debug) {
        MITK_INFO << "Saving connected threshold.";
        mitk::IOUtil::Save(s2c, debugDir + "/seg_s2c.nii");
    }

    UpdateStep(s2c, SegmentationStep::S2C);

    return s2c;
}

mitk::Image::Pointer CemrgFourChamberTools::S2D(mitk::Image::Pointer aorta, 
                                                mitk::Image::Pointer PArt, 
                                                mitk::Image::Pointer svc_slicer, 
                                                mitk::Image::Pointer ivc_slicer,
                                                int aortaSlicerLabel,
                                                int PArtSlicerLabel) {
    if (currentStep != SegmentationStep::S2C) {
        MITK_WARN << "Need to run S2C first";
        return nullptr;
    }
    // aortaSlicerLabel,PArtSlicerLabel
    int aortaLabel = chosenLabels.Get(LabelsType::AORTA);
    int PArtLabel = chosenLabels.Get(LabelsType::PULMONARY_ARTERY);
    int svcLabel = chosenLabels.Get(LabelsType::SVC);
    int ivcLabel = chosenLabels.Get(LabelsType::IVC);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(aortaLabel);

    mitk::Image::Pointer s2d = AddMaskReplaceOnly(currentImage, aorta, aortaSlicerLabel, labelsToProcess);

    labelsToProcess.clear();
    labelsToProcess.push_back(PArtLabel);

    s2d = AddMaskReplaceOnly(s2d, PArt, PArtSlicerLabel, labelsToProcess);

    mitk::Image::Pointer newRA = LabelMaskAndOperation(s2d, svc_slicer, svcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(svcLabel);
    s2d = AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    newRA = LabelMaskAndOperation(s2d, ivc_slicer, ivcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(ivcLabel);
    s2d = AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    if (debug) { 
        MITK_INFO << "Saving s2d.";
        mitk::IOUtil::Save(s2d, debugDir + "/seg_s2d.nii");
    }

    return s2d;
    UpdateStep(s2d, SegmentationStep::S2D);

}

mitk::Image::Pointer CemrgFourChamberTools::S2E(std::vector<unsigned int> seedSVC){
    if (currentStep != SegmentationStep::S2D) {
        MITK_WARN << "Need to run S2D first";
        return nullptr;
    }

    int svcLabel = chosenLabels.Get(LabelsType::SVC);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    mitk::Image::Pointer s2e = ConnectedComponent(currentImage, seedSVC, svcLabel);

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(svcLabel);
    s2e = AddMaskReplaceOnly(s2e, s2e, RALabel, labelsToProcess);

    mitk::Image::Pointer cc = mitk::IOUtil::Load<mitk::Image>(debugDir + "/CC.nii");
    if (!cc) {
        MITK_WARN << "Could not load CC.nii";
        return nullptr;
    }

    s2e = AddMaskReplace(s2e, cc, svcLabel);

    if (debug) { 
        MITK_INFO << "Saving s2e.";
        mitk::IOUtil::Save(s2e, debugDir + "/seg_s2e.nii");
    }

    UpdateStep(s2e, SegmentationStep::S2E);
    return s2e;
}

mitk::Image::Pointer CemrgFourChamberTools::S2F(std::vector<unsigned int> seedIVC) {
    if (currentStep != SegmentationStep::S2E) {
        MITK_WARN << "Need to run S2E first";
        return nullptr;
    }

    int ivcLabel = chosenLabels.Get(LabelsType::IVC);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    mitk::Image::Pointer s2f = ConnectedComponent(currentImage, seedIVC, ivcLabel);

    mitk::Image::Pointer cc = mitk::IOUtil::Load<mitk::Image>(debugDir + "/CC.nii");
    if (!cc) {
        MITK_WARN << "Could not load CC.nii";
        return nullptr;
    }

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(ivcLabel);
    s2f = AddMaskReplaceOnly(s2f, s2f, RALabel, labelsToProcess);
    s2f = AddMaskReplace(s2f, cc, ivcLabel);

    if (debug) { 
        MITK_INFO << "Saving s2f.";
        mitk::IOUtil::Save(s2f, debugDir + "/seg_s2f.nii");
    }

    UpdateStep(s2f, SegmentationStep::S2F);

    return s2f;
}


mitk::Image::Pointer CemrgFourChamberTools::CropSvcIvc(std::vector<mitk::Image::Pointer> images,
                                                       std::vector<unsigned int> seedSVC,
                                                       std::vector<unsigned int> seedIVC, 
                                                       int aortaSlicerLabel,  
                                                       int PArtSlicerLabel) {
    if (images.size() != 5) {
        MITK_WARN << "Need 5 images to crop SVC/IVC: seg_s2a, aorta_slicer, PArt_slicer, SVC_slicer and IVC_slicer images";
        return nullptr;
    }

    int svcLabel = chosenLabels.Get(LabelsType::SVC);
    int ivcLabel = chosenLabels.Get(LabelsType::IVC);
    int aortaLabel = chosenLabels.Get(LabelsType::AORTA);
    int PArtLabel = chosenLabels.Get(LabelsType::PULMONARY_ARTERY);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    // s2a = images.at(0);
    mitk::Image::Pointer s2b = ConnectedComponentKeep(images.at(0), seedSVC, svcLabel);
    if (debug) {
        MITK_INFO << "Saving connected threshold.";
        mitk::IOUtil::Save(s2b, debugDir + "/seg_s2b.nii");
    }

    mitk::Image::Pointer s2c = ConnectedComponentKeep(s2b, seedIVC, ivcLabel);
    if (debug) {
        MITK_INFO << "Saving connected threshold.";
        mitk::IOUtil::Save(s2c, debugDir + "/seg_s2c.nii");
    }

    std::vector<int> labelsToProcess; 
    labelsToProcess.push_back(aortaLabel);

    // aorta = images.at(1);
    mitk::Image::Pointer s2d = AddMaskReplaceOnly(s2c, images.at(1), aortaSlicerLabel, labelsToProcess);

    labelsToProcess.clear();
    labelsToProcess.push_back(PArtLabel);

    // PArt = images.at(2);
    s2d = AddMaskReplaceOnly(s2d, images.at(2), PArtSlicerLabel, labelsToProcess);

    // svc = images.at(3);
    mitk::Image::Pointer newRA = LabelMaskAndOperation(s2d, images.at(3), svcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(svcLabel);
    s2d = AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    
    // ivc = images.at(4);
    newRA = LabelMaskAndOperation(s2d, images.at(4), ivcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(ivcLabel);
    s2d = AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    if (debug) { 
        MITK_INFO << "Saving s2d.";
        mitk::IOUtil::Save(s2d, debugDir + "/seg_s2d.nii");
    }

    mitk::Image::Pointer s2e = ConnectedComponent(s2d, seedSVC, svcLabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(svcLabel);
    s2e = AddMaskReplaceOnly(s2e, s2e, RALabel, labelsToProcess);

    mitk::Image::Pointer cc = mitk::IOUtil::Load<mitk::Image>(debugDir + "/CC.nii");
    s2e = AddMaskReplace(s2e, cc, svcLabel);

    if (debug) { 
        MITK_INFO << "Saving s2e.";
        mitk::IOUtil::Save(s2e, debugDir + "/seg_s2e.nii");
    }

    mitk::Image::Pointer s2f = ConnectedComponent(s2e, seedIVC, ivcLabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(ivcLabel);
    s2f = AddMaskReplaceOnly(s2f, s2f, RALabel, labelsToProcess);
    s2f = AddMaskReplace(s2f, cc, ivcLabel);

    return s2f;
}

mitk::Image::Pointer CemrgFourChamberTools::ConnectedComponent(mitk::Image::Pointer seg, std::vector<unsigned int> seedIdx, int layer, bool keep) {
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
    mitk::IOUtil::Save(ccImage, debugDir + "/CC.nii");

    int replaceValue = 0;
    if (keep) {
        seg = RemoveLabel(seg, layer);
        replaceValue = layer;
    }

    return AddMaskReplace(seg, ccImage, replaceValue);

}
