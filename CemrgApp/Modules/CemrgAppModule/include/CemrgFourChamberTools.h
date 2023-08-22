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

#include "CemrgMultilabelSegmentationUtils.h"
#include "FourChamberCommon.h"
#include "CemrgCommandLine.h"

struct Step {
    std::string name;
    std::string imageFilename;
    mitk::Image::Pointer image;

    Step(std::string stepName) : name(stepName), imageFilename("") {
        image = nullptr;
    }

    inline void SetImageFileName(std::string dir) { imageFilename = dir + "/" + name + ".nii"; };
    inline std::string GetImageFileName() { return imageFilename; };

    inline void SetImage(mitk::Image::Pointer img) { image = img; };
    inline void DropImage() { image = nullptr; }; 
    inline void SaveImage() { mitk::IOUtil::Save(image, imageFilename); };
    inline void LoadImage() { image = mitk::IOUtil::Load<mitk::Image>(imageFilename); };
};

struct Stage {
    std::string name;
    std::vector<Step> steps;

    Stage() : name(""){};
    Stage(std::string stageName) : name(stageName){};

    inline int StepsInStage() { return steps.size(); };
};

class SegmentationStepManager {
public:
    SegmentationStepManager();
    SegmentationStepManager(std::string path);
    ~SegmentationStepManager();

    void SetPathToImages(std::string path);

    void NextStep();
    void NextStage();
    inline Step* GetCurrentStep() { return &(stages[stageIndx].steps[stepIndx]); };
    inline Stage* GetCurrentStage() { return &(stages[stageIndx]); };
    void SaveCurrentStepImage();
    std::string PrintCurrentStepInfo();

    inline QString GetStepName() { return QString::fromStdString( GetCurrentStep()->name); };
    inline mitk::Image::Pointer GetStepImage() { return GetCurrentStep()->image; };
    inline void LoadImageForCurrentStep() { GetCurrentStep()->LoadImage(); };

    void UpdateStepWithImage(mitk::Image::Pointer image, bool debug);

    void NavigateToStep(const std::string &stageName, const std::string &stepName);
    void NavigateToStepFromFile(const QString& filename);
    QJsonObject GetJson();
    void SaveToJson(const QString& filename);
    void ModifyJson(); 

private:
    std::vector<Stage> stages;
    std::string pathToImages, pathStepsJson;
    
    size_t stageIndx;
    size_t stepIndx;
    bool started;
};

class MITKCEMRGAPPMODULE_EXPORT CemrgFourChamberTools {
    typedef itk::Image<uint8_t, 3> ImageType;
    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    
    public:
        CemrgFourChamberTools();
        ~CemrgFourChamberTools();

        inline void SetDebug(bool value) { debug = value; };
        inline void SetDebugOn() { SetDebug(true); };
        inline void SetDebugOff() { SetDebug(false); };

        void SetSegDir(std::string segDirectory);

        inline std::string GetSegDir() { return segDir; };
        inline std::string GetDebugDir() { return directory; };
        inline QString QGetSegDir() { return QString::fromStdString(segDir); };
        inline QString QGetDebugDir() { return QString::fromStdString(directory); };

        inline void SetLabel(LabelsType stt, int label) { chosenLabels.Set(stt, label); };

        // tools to manipulate segStepManager
        void UpdateSegmentationStep(mitk::Image::Pointer newImage);
        void SaveSegmentationStage(QString pathname) { segStepManager.SaveToJson(pathname); };
        inline std::string StepName() { return segStepManager.GetStepName().toStdString(); };
        inline std::string StepFileImage() { return segStepManager.GetCurrentStep()->GetImageFileName(); };
        inline void UpdateCurrentImage() { currentImage = segStepManager.GetStepImage(); };
        inline void NavigateToStep(QString filepath) { segStepManager.NavigateToStepFromFile(filepath); };
        inline mitk::Image::Pointer GetCurrentImage() { return currentImage; };

        // Manipulation utilities
        inline void SetCylinders(QJsonObject json){ _cylinders.SetPointsFromJson(json); };
        inline void SetSlicers(QJsonObject json) { _slicers.SetPointsFromJson(json); };
        inline void SetValvePoints(QJsonObject json) { _valvePoints.SetPointsFromJson(json); };

        inline bool CylindersSet() { return _cylinders.IsPointSet(); };
        inline bool SlicersSet() { return _slicers.IsPointSet(); };
        inline bool ValvePointsSet() { return _valvePoints.IsPointSet(); };


        mitk::Image::Pointer Cylinder(mitk::Image::Pointer seg, QString ptPrefix, double slicerRadius, double slicerHeight, ManualPoints mpl, QString saveAs = "");
        mitk::Image::Pointer CreateSvcIvc(std::vector<mitk::Image::Pointer> images, int RspvLabel=10, int SvcLabel=13, int IvcLabel=14);
        mitk::Image::Pointer CropSvcIvc(std::vector<mitk::Image::Pointer> images, int aortaSlicerLabel = 0 , int PArtSlicerLabel = 0);

    protected:
        mitk::Image::Pointer S2(std::vector<unsigned int> seeds, QString checkPrevious, LabelsType lt);
        inline mitk::Image::Pointer S2B(std::vector<unsigned int> seedSVC) { return S2(seedSVC, "s2a", LabelsType::SVC); };
        inline mitk::Image::Pointer S2C(std::vector<unsigned int> seedIVC) { return S2(seedIVC, "s2b", LabelsType::IVC); };

        mitk::Image::Pointer S2D(mitk::Image::Pointer aorta, mitk::Image::Pointer PArt, mitk::Image::Pointer svc_slicer, mitk::Image::Pointer ivc_slicer, int aortaSlicerLabel, int PArtSlicerLabel);
        mitk::Image::Pointer S2E(std::vector<unsigned int> seedSVC);
        mitk::Image::Pointer S2F(std::vector<unsigned int> seedIVC);

    private : 
        CylinderPointsType _cylinders;
        SlicersPointsType _slicers;
        ValvePlainsPointsType _valvePoints;

        bool debug;
        std::string segDir, directory;

        mitk::Image::Pointer currentImage;
        SegmentationStepManager segStepManager;
        LabelsStruct chosenLabels;
        SegmentationLabels segmentationLabels;
};
#endif // CemrgFourChamberTools_h
