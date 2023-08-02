/*=========================================================================
BSD 3-Clause License

Copyright (c) 2020, CemrgAppDevelopers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/


#ifndef FourChamberView_h
#define FourChamberView_h

#include <berryISelectionListener.h>

#include <QmitkAbstractView.h>

#include "ui_FourChamberViewControls.h"
#include "ui_Meshtools3DParameterUI.h"

#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

#include "CemrgFourChamberTools.h"
#include "FourChamberCommon.h"
#include "CemrgDataInteractor.h"

/**
  \brief FourChamberView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class FourChamberView : public QmitkAbstractView {
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:
    static const std::string VIEW_ID;
    static const QString POINTS_FILE;   // all the user-selected points
    static const QString POINTS_FILE_INDEX;   // all the user-selected points
    static const QString GEOMETRY_FILE; // origin and spacing

    // helper functions
    bool RequestAnyFolderFromUser(QString & dir, std::string msg, bool project_dir=false);
    bool CheckForExistingFile(QString dir, QString filename);
    int Ask(std::string title, std::string msg);
    void Warn(std::string title, std::string msg);
    void Inform(std::string title, std::string msg);
    void PrintMeshingParameters(QString path_to_par);
    void LoadMeshingParametersFromJson(QString dir, QString json_file);
    QString GetPointTypeString(ManualPointsType mpt);
    QString ArrayToString(double *arr, int size, QString title);
    bool ArrayEqual(double *arr1, double *arr2, int size, double tol=0.0001);
    void ParseJsonArray(QJsonObject json, QString key, double *arr, int size=3);
    bool CheckPointsInJsonFile(ManualPointsType mpt);
    void ReloadJsonPoints();
    std::string PrintPoints(QJsonObject json, QStringList keysList, QString title);

    void
    SetButtonsEnable(bool enable);
    inline void SetButtonsEnableToOn(){ SetButtonsEnable(true); };
    inline void SetButtonsEnableToOff() { SetButtonsEnable(false); };

    void InitialiseJsonObjects();
    QStringList GetPointLabelOptions(ManualPointsType mpt);
    CemrgDataInteractor::Pointer CreateInteractorWithOptions(ManualPointsType mpt);
    void InitialiseQStringListsFromSize(int num, QStringList &values, QStringList &types);

    // User Select Functions
    bool UserSelectMeshtools3DParameters(QString pre_input_path);

    // inline means they're defined here, not in the cpp filemguvc
    inline bool RequestProjectDirectoryFromUser(){ return RequestAnyFolderFromUser(directory, "Project folder", true); };
    inline bool RequestCarpDirectoryFromUser(){ return RequestAnyFolderFromUser(carp_directory, "CARP binary folder", false); };
    inline QString Path(QString fnameExt = "") {return (directory + "/" + fnameExt); };
    inline std::string StdStringPath(QString fnameExt=""){return (Path(fnameExt).toStdString());};
    

protected slots:
    // here you add the functions that willl be linked to buttons
    void SetWorkingFolder();
    void PrepareSegmentation();
    void Meshing();
    void UVCs();
    void VentricularFibres();
    void AtrialFibres();
    void SimulationSetup();

    void LoadDICOM();
    void GetOriginSpacing();
    void SegmentImgs();
    void ExtractSurfaces();
    void SelectLARALandmarks();
    void CalculateUVCs();

    void Corrections();
    void SelectPoints();
    void SelectPointsCylinders();
    void SelectPointsSlicers();
    void SelectPointsValvePlains();
    void SelectPointsCheck();
    void SelectPointsReset();

    void M3dBrowseFile(const QString &dir);

protected:
    // this whole block hardly ever changes 
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::FourChamberViewControls m_Controls;
    Ui::Meshtools3DParameterUI m_m3d;

private:
    // put here the things which belong to the class, like working folder name, etc
    QString fileName, directory, current_seg_name, carp_directory;
    QStringList pt_keys_init, pt_keys_slicers, pt_keys_final;
    QJsonObject json_points, json_geometry; // keeps all the points available
    bool points_file_loaded; // keeps track if points.json has been loaded
    bool carpless;           // true if user does not have CARP installed

    mitk::PointSet::Pointer pset;
    std::unique_ptr<CemrgFourChamberTools> fourch_object;

    M3DParameters meshing_parameters;
    FourChamberSubfolders SDIR; // helps set access subfolders in working directory
    SegmentationPointsIds SegPointIds; // keeps track of the segmentation points
    double origin[3], spacing[3];
};

#endif // FourChamberView_h
