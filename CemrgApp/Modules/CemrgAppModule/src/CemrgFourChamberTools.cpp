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
#include <QJsonArray>

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

SegmentationStepManager::SegmentationStepManager() : stageIndx(0), stepIndx(0), started(false), pathStepsJson("") {
    // Define and populate stages and steps
    Stage segmentationStage("Segmentation");
    segmentationStage.steps.push_back(Step("s1"));

    Stage cylindersStage("Cylinders");
    QStringList stgs = QStringList() << "a" << "b"  << "c" << "d" << "e" << "f";
    for (int ix = 0; ix < stgs.size(); ix++) {
       cylindersStage.steps.push_back(Step(("s2" + stgs.at(ix)).toStdString()));
    }
    stages.push_back(cylindersStage);

    Stage slicersStage("Slicers");
    stgs << "g" << "h"  << "i" << "j" << "k" << "l" << "m" << "n" << "o" << "p";
    for (int ix = 0; ix < stgs.size(); ix++) {
       slicersStage.steps.push_back(Step(("s3" + stgs.at(ix)).toStdString()));
    }
    stages.push_back(slicersStage);

    stgs.clear();

    Stage valvePlainsStage("ValvePlains");
    stgs << "s3r" << "s3s" << "s4a" << "s4b" << "s4c" << "s4d" 
        << "s4e" << "s4f" << "s4g" << "s4h" << "s4i" << "s4j" << "s4k";
    for (int ix = 0; ix < stgs.size(); ix++) {
       valvePlainsStage.steps.push_back(Step(stgs.at(ix).toStdString()));
    }
    stages.push_back(valvePlainsStage);

    Stage cleanSegmentationStage("CleanSegmentation");
    cleanSegmentationStage.steps.push_back(Step("s5"));

    stages.push_back(cleanSegmentationStage);
}

SegmentationStepManager::SegmentationStepManager(std::string path) {
    SegmentationStepManager();
    SetPathToImages(path);
}

SegmentationStepManager::~SegmentationStepManager() {
    for (size_t ix = 0; ix < stages.size(); ix++) {
        for (size_t jx = 0; jx < stages.at(ix).steps.size(); jx++) {
            stages.at(ix).steps.at(jx).DropImage();
        }
    }
}

void SegmentationStepManager::SetPathToImages(std::string path) {
    for (size_t ix = 0; ix < stages.size(); ix++) {
        for (size_t jx = 0; jx < stages.at(ix).steps.size(); jx++) {
            stages.at(ix).steps.at(jx).SetImageFileName(path);
        }
    }
}

void SegmentationStepManager::NextStep()
{
    if (stages[stageIndx].steps[stepIndx].image) {
        stages[stageIndx].steps[stepIndx].DropImage();
    }

    stepIndx++;
    if (stepIndx >= stages[stageIndx].steps.size()) {
        NextStage();
        stepIndx = 0;
    }
}

void SegmentationStepManager::NextStage() {
    stageIndx++;
    MITK_INFO(stageIndx >= stages.size()) << "Segmentation steps complete";
}

void SegmentationStepManager::SaveCurrentStepImage() {
    GetCurrentStep()->SaveImage();
}

std::string SegmentationStepManager::PrintCurrentStepInfo() {
    std::string msg = "Stage: " + stages[stageIndx].name + "\n";
    msg += "Step: " + stages[stageIndx].steps[stepIndx].name;
    return msg;
}

void SegmentationStepManager::UpdateStepWithImage(mitk::Image::Pointer image, bool debug) {
    if (!image) {
        MITK_WARN << "Image for update was NULL. Step not updated.";
        return;
    }

    if (started) {
        NextStep();
    } else {
        started = true;
    }
    
    GetCurrentStep()->SetImage(image);

    if (debug) {
        SaveCurrentStepImage();
    }

}

void SegmentationStepManager::NavigateToStep(const std::string &stageName, const std::string &stepName) {
    for (size_t stagex = 0; stagex < stages.size(); ++stagex) {
        if (stages.at(stagex).name == stageName) {
            stageIndx = stagex;
            for (size_t stepx = 0; stepx < stages.at(stagex).steps.size(); ++stepx) {
                if (stages.at(stagex).steps.at(stepx).name == stepName) {
                    stepIndx = stepx;
                    return;
                }
            }
        }
        return;
    }
}

void SegmentationStepManager::NavigateToStepFromFile(QString filename) {
    
    QFileInfo fi(filename);
    QJsonObject json = CemrgCommonUtils::ReadJSONFile(fi.absolutePath(), fi.fileName());

    pathToImages = json["pathToImages"].toString().toStdString();
    started = json["started"].toBool();
    stageIndx = json["stageIndx"].toInt();
    stepIndx = json["stepIndx"].toInt();

    GetCurrentStep()->SetImageFileName(pathToImages);
    GetCurrentStep()->LoadImage();

    MITK_INFO << PrintCurrentStepInfo();
}

void SegmentationStepManager::SaveToJson(QString filename) {
    QJsonObject json = GetJson();
    QFileInfo fi(filename);
    CemrgCommonUtils::WriteJSONFile(json, fi.absolutePath(), fi.baseName()+".json");
    pathStepsJson = filename.toStdString();
}

QJsonObject SegmentationStepManager::GetJson() {
    QJsonObject json;
    QJsonArray stagesArray;
    for (const Stage& stage : stages) {
        QJsonObject stageJson;
        stageJson["name"] = QString::fromStdString(stage.name);

        QJsonArray stepsArray;
        for (const Step& step : stage.steps) {
            QJsonObject stepJson;
            stepJson["name"] = QString::fromStdString(step.name);
            stepJson["imageFilename"] = QString::fromStdString(step.imageFilename);
            stepsArray.append(stepJson);
        }
        stageJson["steps"] = stepsArray;
        stagesArray.append(stageJson);
    }

    json["stages"] = stagesArray;
    json["pathToImages"] = QString::fromStdString(pathToImages);
    json["stageIndx"] = QJsonValue::fromVariant( (int) stageIndx);
    json["stepIndx"] = QJsonValue::fromVariant( (int) stepIndx);

    QString boolstr = (started) ? "true" : "false";
    json["started"] = QJsonValue::fromVariant(boolstr);

    return json;
}

void SegmentationStepManager::ModifyJson() {
    QJsonObject json = GetJson();
    SaveToJson(QString::fromStdString(pathStepsJson));
}

// CemrgFourChamberTools
CemrgFourChamberTools::CemrgFourChamberTools(){
    _cylinders = CylinderPointsType();
    _slicers = SlicersPointsType();
    _valvePoints = ValvePlainsPointsType();
    debug = false;
    directory = "";
    segStepManager = SegmentationStepManager();
}

CemrgFourChamberTools::~CemrgFourChamberTools(){
    segStepManager.~SegmentationStepManager();
}

void CemrgFourChamberTools::SetSegDir(std::string segDirectory) {
    segDir = segDirectory;
    directory = segDir + "/tmp";
    segStepManager.SetPathToImages(directory);  
}

void CemrgFourChamberTools::UpdateSegmentationStep(mitk::Image::Pointer newImage) {
    MITK_INFO << "Updating segmentation step: " + segStepManager.PrintCurrentStepInfo();
    segStepManager.UpdateStepWithImage(newImage, debug);
    UpdateCurrentImage();
}

void CemrgFourChamberTools::SaveSegmentationStage(QString pathname) { 
    segStepManager.SaveToJson(pathname); 
}

void CemrgFourChamberTools::NavigateToStep(QString filepath) { 
    segStepManager.NavigateToStepFromFile(filepath); 
}

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
    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    mitk::Image::Pointer svcIvc = segUtils->AddMaskReplaceOnly(seg, svc, SvcLabel, labelsToProcess);
    if (svcIvc == nullptr) return nullptr;

    svcIvc = segUtils->AddMask(svcIvc, ivc, IvcLabel);
    if (svcIvc == nullptr) return nullptr;

    return svcIvc;
}

mitk::Image::Pointer CemrgFourChamberTools::S2(std::vector<unsigned int> seeds, QString checkPrevious, LabelsType lt){
    // used for S2B and S2C
    if (segStepManager.GetStepName() != checkPrevious) {
        MITK_WARN << ("Need to run " + checkPrevious.toUpper() + " first").toStdString();
        return nullptr;
    }

    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    
    mitk::Image::Pointer s2 = segUtils->ConnectedComponent(segStepManager.GetStepImage(), seeds, chosenLabels.Get(lt));
    segStepManager.UpdateStepWithImage(s2, debug);

    return s2;

}

mitk::Image::Pointer CemrgFourChamberTools::S2D(mitk::Image::Pointer aorta,  mitk::Image::Pointer PArt,  mitk::Image::Pointer svc_slicer,  mitk::Image::Pointer ivc_slicer, int aortaSlicerLabel, int PArtSlicerLabel) {
    if (segStepManager.GetStepName() != "s2c") {
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

    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    mitk::Image::Pointer s2d = segUtils->AddMaskReplaceOnly(segStepManager.GetStepImage(), aorta, aortaSlicerLabel, labelsToProcess);

    labelsToProcess.clear();
    labelsToProcess.push_back(PArtLabel);

    s2d = segUtils->AddMaskReplaceOnly(s2d, PArt, PArtSlicerLabel, labelsToProcess);

    mitk::Image::Pointer newRA = segUtils->LabelMaskAndOperation(s2d, svc_slicer, svcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(svcLabel);
    s2d = segUtils->AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    newRA = segUtils->LabelMaskAndOperation(s2d, ivc_slicer, ivcLabel, RALabel);

    labelsToProcess.clear();
    labelsToProcess.push_back(ivcLabel);
    s2d = segUtils->AddMaskReplaceOnly(s2d, newRA, RALabel, labelsToProcess);

    segStepManager.UpdateStepWithImage(s2d, debug);

    return s2d;

}

mitk::Image::Pointer CemrgFourChamberTools::S2E(std::vector<unsigned int> seedSVC){
    if (segStepManager.GetStepName() != "s2d") {
        MITK_WARN << "Need to run S2D first";
        return nullptr;
    }

    int svcLabel = chosenLabels.Get(LabelsType::SVC);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    segUtils->SetDebugDir(QString::fromStdString(directory));
    mitk::Image::Pointer s2e = segUtils->ConnectedComponent(segStepManager.GetStepImage(), seedSVC, svcLabel);

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(svcLabel);
    s2e = segUtils->AddMaskReplaceOnly(s2e, s2e, RALabel, labelsToProcess);

    mitk::Image::Pointer cc = mitk::IOUtil::Load<mitk::Image>(directory + "/CC.nii");
    if (!cc) {
        MITK_WARN << "Could not load CC.nii";
        return nullptr;
    }

    s2e = segUtils->AddMaskReplace(s2e, cc, svcLabel);
    segStepManager.UpdateStepWithImage(s2e, debug);

    return s2e;
}

mitk::Image::Pointer CemrgFourChamberTools::S2F(std::vector<unsigned int> seedIVC) {
    if (segStepManager.GetStepName() != "s2e") {
        MITK_WARN << "Need to run S2E first";
        return nullptr;
    }

    int ivcLabel = chosenLabels.Get(LabelsType::IVC);
    int RALabel = chosenLabels.Get(LabelsType::RIGHT_ATRIUM);

    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    mitk::Image::Pointer s2f = segUtils->ConnectedComponent(segStepManager.GetStepImage(), seedIVC, ivcLabel);

    mitk::Image::Pointer cc = mitk::IOUtil::Load<mitk::Image>(directory + "/CC.nii");
    if (!cc) {
        MITK_WARN << "Could not load CC.nii";
        return nullptr;
    }

    std::vector<int> labelsToProcess;
    labelsToProcess.push_back(ivcLabel);
    s2f = segUtils->AddMaskReplaceOnly(s2f, s2f, RALabel, labelsToProcess);
    s2f = segUtils->AddMaskReplace(s2f, cc, ivcLabel);

    segStepManager.UpdateStepWithImage(s2f, debug);

    return s2f;
}

mitk::Image::Pointer CemrgFourChamberTools::CropSvcIvc(std::vector<mitk::Image::Pointer> images, int aortaSlicerLabel, int PArtSlicerLabel) {

    std::vector<double> seedSVCWorld(3), seedIVCWorld(3);
    std::vector<unsigned int> seedSVC(3), seedIVC(3);

    for (int ix = 0; ix < 3; ix++) {
        seedSVCWorld.at(ix) = _slicers.GetPointAt(PointsNamesType::SVC_TIP, ix);
        seedIVCWorld.at(ix) = _slicers.GetPointAt(PointsNamesType::IVC_TIP, ix);
    }

    std::unique_ptr<CemrgMultilabelSegmentationUtils> segUtils = std::unique_ptr<CemrgMultilabelSegmentationUtils>(new CemrgMultilabelSegmentationUtils());
    segUtils->WorldToIndex(images.at(0), seedSVCWorld, seedSVC);
    segUtils->WorldToIndex(images.at(0), seedIVCWorld, seedIVC);


    if (images.size() != 5) {
        MITK_WARN << "Need 5 images to crop SVC/IVC: seg_s2a, aorta_slicer, PArt_slicer, SVC_slicer and IVC_slicer images";
        return nullptr;
    }

    if (!S2B(seedSVC)) {
        MITK_WARN << "Could not run S2B";
        return nullptr;
    }

    if (!S2C(seedIVC)) {
        MITK_WARN << "Could not run S2C";
        return nullptr;
    }

    if (!S2D(images.at(1), images.at(2), images.at(3), images.at(4), aortaSlicerLabel, PArtSlicerLabel)) {
        MITK_WARN << "Could not run S2D";
        return nullptr;
    }

    if (!S2E(seedSVC)){
        MITK_WARN << "Could not run S2E";
        return nullptr;
    }

    return S2F(seedIVC);
}