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

#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

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

    // helper functions
    bool RequestProjectDirectoryFromUser();
    int Ask(std::string title, std::string msg);

    void InitialiseJsonObjects();
    QStringList GetPointLabelOptions(QString opt);
    void CreateInteractorWithOptions(QString opts);
    void InitialiseQStringListsFromSize(int num, QStringList &values, QStringList &types);

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

    void SelectPointsInit();
    void SelectPointsSlicers();
    void SelectPointsFinal();

protected:
    // this whole block hardly ever changes 
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::FourChamberViewControls m_Controls;

private:
    // put here the things which belong to the class, like working folder name, etc
    QString fileName, directory, current_seg_name;
    QStringList pt_keys_init, pt_keys_slicers, pt_keys_final;
    QJsonObject json_init, json_slicers, json_final;
};

#endif // FourChamberView_h
