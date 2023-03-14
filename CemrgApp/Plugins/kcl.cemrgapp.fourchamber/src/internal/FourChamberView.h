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

/// @brief Meshtool3D static libraries parameters structures
struct M3DParameters {
    QString seg_dir;
    QString seg_name;
    bool mesh_from_segmentation;
    bool boundary_relabelling;

    double facet_angle;
    double facet_size;
    double facet_distance;
    double cell_rad_edge_ratio;
    double cell_size;
    double rescaleFactor;

    double abs_tol;
    double rel_tol;
    int itr_max;
    int dimKrilovSp;
    bool verbose;

    bool eval_thickness;

    QString out_dir;
    QString out_name;

    bool out_medit;
    bool out_carp;
    bool out_carp_binary;
    bool out_vtk;
    bool out_vtk_binary;
    bool out_potential;

    M3DParameters()
    :seg_dir(""),
     seg_name(""),
     mesh_from_segmentation(true),
     boundary_relabelling(false),
     facet_angle(30.0),
     facet_size(0.8),
     facet_distance(4),
     cell_rad_edge_ratio(2.0),
     cell_size(0.8),
     rescaleFactor(1000),
     abs_tol(1e-6),
     rel_tol(1e-6),
     itr_max(700),
     dimKrilovSp(500),
     verbose(true),
     eval_thickness(false),
     out_dir(""),
     out_name(""),
     out_medit(false),
     out_carp(true),
     out_carp_binary(false),
     out_vtk(false),
     out_vtk_binary(false),
     out_potential(false)
    {}

    bool CheckFormats(){
        return !(out_medit || out_carp || out_carp_binary ||out_vtk || out_vtk_binary);
    }

    void SetBinaryFormats(){
        out_carp_binary = out_carp;
        out_vtk_binary = out_vtk;
    }

    void KeysAndValues(QStringList& keys, QStringList& values, QStringList& types){
        keys << "seg_dir" << "seg_name" << "mesh_from_segmentation" << "boundary_relabelling" << 
                "facet_angle" << "facet_size" << "facet_distance" << "cell_rad_edge_ratio" << "cell_size" << "rescaleFactor" << 
                "abs_tol" << "rel_tol" << "itr_max" << "dimKrilovSp" << "verbose" << "eval_thickness" << 
                "out_dir" << "out_name" << "out_medit" << "out_carp" << "out_carp_binary" << "out_vtk" << "out_vtk_binary" << "out_potential";
        values << seg_dir << seg_name << QString::number(mesh_from_segmentation) << QString::number(boundary_relabelling)
               << QString::number(facet_angle) << QString::number(facet_size) << QString::number(facet_distance) << QString::number(cell_rad_edge_ratio) 
               << QString::number(cell_size) << QString::number(rescaleFactor) << QString::number(abs_tol)
               << QString::number(rel_tol) << QString::number(itr_max) << QString::number(dimKrilovSp) << QString::number(verbose) << QString::number(eval_thickness) 
               << out_dir << out_name
               << QString::number(out_medit) << QString::number(out_carp) << QString::number(out_carp_binary) << QString::number(out_vtk) << QString::number(out_vtk_binary) << QString::number(out_potential);
        types << "string" << "string" << "int" << "int" << 
                "float" << "float" << "float" << "float" << "float" << "int" << 
                "float" << "float" << "int" << "int" << "int" << "int" << 
                "string" << "string" << "int" << "int" << "int" << "int" << "int" << "int";
    }

};


enum ManualPointsType { CYLINDERS, SLICERS, VALVE_PLAINS };
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
    static const QString POINTS_FILE;
    static const QString MESH_SDIR;
    static const QString SEG_SDIR;
    // helper functions
    bool RequestProjectDirectoryFromUser();
    int Ask(std::string title, std::string msg);
    void Warn(std::string title, std::string msg);
    void PrintMeshingParameters(QString path_to_par);
    void LoadMeshingParametersFromJson(QString dir, QString json_file);
    QString GetPointTypeString(ManualPointsType mpt);

    void InitialiseJsonObjects();
    QStringList GetPointLabelOptions(ManualPointsType mpt);
    void CreateInteractorWithOptions(ManualPointsType mpt);
    void InitialiseQStringListsFromSize(int num, QStringList &values, QStringList &types);

    // User Select Functions
    bool UserSelectMeshtools3DParameters(QString pre_input_path);

    // inline means they're defined here, not in the cpp file
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
    void ExtractSurfaces();
    void SelectLARALandmarks();
    void CalculateUVCs();

    void SelectPointsCylinders();
    void SelectPointsSlicers();
    void SelectPointsValvePlains();

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
    QString fileName, directory, current_seg_name;
    QStringList pt_keys_init, pt_keys_slicers, pt_keys_final;
    QJsonObject json_points; // keeps all the points available
    bool points_file_loaded; // keeps track if points.json has been loaded

    M3DParameters meshing_parameters;
};

#endif // FourChamberView_h
