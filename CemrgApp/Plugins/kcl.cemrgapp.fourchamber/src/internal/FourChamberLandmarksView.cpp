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
 * Morphological Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include "FourChamberLandmarksView.h"
#include "FourChamberView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkClipPolyData.h>

//ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

// Qt
#include <QMessageBox>
#include <QDesktopWidget>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>

QString FourChamberLandmarksView::directory;
QString FourChamberLandmarksView::laSubdir;
QString FourChamberLandmarksView::raSubdir;
QString FourChamberLandmarksView::laName;
QString FourChamberLandmarksView::raName;

const std::string FourChamberLandmarksView::VIEW_ID = "org.mitk.views.fourchamberlandmarksview";

void FourChamberLandmarksView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_help, SIGNAL(clicked()), this, SLOT(Help()));
    connect(m_Controls.button_save, SIGNAL(clicked()), this, SLOT(SaveSelectedPoints()));
    connect(m_Controls.combo_select, SIGNAL(currentIndexChanged(int)), this, SLOT(ComboSelectWorkingMesh(int)));

    //Create GUI widgets
    inputsPickedPoints = new QDialog(0,0);
    m_PickedPoints.setupUi(inputsPickedPoints);
    connect(m_PickedPoints.buttonBox, SIGNAL(accepted()), inputsPickedPoints, SLOT(accept()));
    connect(m_PickedPoints.buttonBox, SIGNAL(rejected()), inputsPickedPoints, SLOT(reject()));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.5,0.5,0.5);
    renderer->AutomaticLightCreationOn();
    renderer->LightFollowCameraOn();
    // renderer->TwoSidedLightingOn();
    // renderer->UpdateLightsGeometryToFollowCamera();
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    std::string shortcuts = "SPACE: Pick point\nDELETE: Remove last point\nH: Help\n";
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->SetRenderWindow(renderWindow);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    callBack->SetCallback(KeyCallBackFunc);
    callBack->SetClientData(this);
    interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    //Initialisation
    iniPreSurf();
}

void FourChamberLandmarksView::SetFocus() {
    m_Controls.button_help->setFocus();
}

void FourChamberLandmarksView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void FourChamberLandmarksView::SetDirectoryFile(const QString directory, const QString laSubdir, const QString laName, const QString raSubdir, const QString raName) {
    FourChamberLandmarksView::directory = directory;
    FourChamberLandmarksView::laSubdir = laSubdir;
    FourChamberLandmarksView::laName = laName;
    FourChamberLandmarksView::raSubdir = raSubdir;
    FourChamberLandmarksView::raName = raName;
}

FourChamberLandmarksView::~FourChamberLandmarksView() {
    inputsPickedPoints->deleteLater();
}

// slots
void FourChamberLandmarksView::Help(){
    std::string msg = "";
    msg += "This view allows the user to pick landmarks on the atrial surface.\n";
    msg += "Choose the atrium from the combo box, then \n";
    msg += "pick points by pressing the SPACE key.";
    QMessageBox::information(NULL, "Help", msg.c_str());
}

void FourChamberLandmarksView::SaveSelectedPoints() {

} 

void FourChamberLandmarksView::ComboSelectWorkingMesh(int index) {
    QString selected = m_Controls.combo_select->itemText(index);
    surface = mitk::IOUtil::Load<mitk::Surface>(directory.toStdString() + "/" + selected.toStdString());    
    if (surface) {
        if (pickedPoint.IsEmpty()) {
            pickedPoint = PickedPointType();
        } else {
            m_PickedPoints.comboBox->clear();
            pickedPoint.Clear();
        }
        
        QStringList names;
        std::vector<int> availableLabels;
        if (selected.contains(laSubdir)) {
            names << "LA_APEX" << "LA_SEPTUM";
            availableLabels = {AtrialLandmarksType::LA_APEX, AtrialLandmarksType::LA_SEPTUM};
        } else if (selected.contains(raSubdir)) {
            names << "RA_APEX" << "RA_SEPTUM" << "RAA_APEX";
            availableLabels = {AtrialLandmarksType::RA_APEX, AtrialLandmarksType::RA_SEPTUM, AtrialLandmarksType::RAA_APEX};
        }
        pickedPoint.SetAvailableLabels(names, availableLabels);
        m_PickedPoints.comboBox->addItems(names);
    }
}


/// slots - end

void FourChamberLandmarksView::iniPreSurf() {    
    //Load the surface
    QString laPath = directory + "/" + laSubdir + "/" + laName;
    QString raPath = directory + "/" + raSubdir + "/" + raName;

    if (QFile::exists(laPath)) {
        m_Controls.combo_select->addItem(laSubdir + "/" + laName);
    }
    if (QFile::exists(raPath)) {
        m_Controls.combo_select->addItem(raSubdir + "/" + raName);
    }

    if (m_Controls.combo_select->count() > 0) {
        m_Controls.combo_select->setEnabled(false);
        QMessageBox::warning(NULL, "Warning", "No surfaces found in the selected directory!");
    } else { 

        Help();
    }
}

void FourChamberLandmarksView::Visualiser(double opacity){
    MITK_INFO << "[Visualiser]";
    double max_scalar = 19;
    double min_scalar = 1;
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    this->maxScalar = max_scalar;
    this->minScalar = min_scalar;

    if (!pickedPoint.IsEmpty()) {
        SphereSourceVisualiser(pickedPoint.GetLineSeeds());
    }

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUseCellData();
    surfMapper->ScalarVisibilityOn();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.5, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(surfMapper->GetLookupTable());
    scalarBar->SetWidth(0.075);
    scalarBar->SetHeight(0.3);
    scalarBar->SetTextPositionToPrecedeScalarBar();
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue( 0.9, 0.01 );
    scalarBar->SetNumberOfLabels(5);

    std::string scalarBarTitle = "Atrium Labels";
    scalarBar->SetTitle(scalarBarTitle.c_str());
    scalarBar->SetVerticalTitleSeparation(15);

    surfMapper->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(opacity);
    renderer->AddActor(surfActor);
    renderer->AddActor2D(scalarBar);
}

void FourChamberLandmarksView::SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour, double scaleFactor){
    // e.g colour = "0.4,0.1,0.0" - values for r,g, and b separated by commas.
    double r, g, b;
    bool okr, okg, okb;
    r = colour.section(',',0,0).toDouble(&okr);
    g = colour.section(',',1,1).toDouble(&okg);
    b = colour.section(',',2,2).toDouble(&okb);

    if(!okr) r=1.0;
    if(!okg) g=0.0;
    if(!okb) b=0.0;

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkSphereSource> glyphSource = vtkSmartPointer<vtkSphereSource>::New();
    glyph3D->SetInputData(pointSources);
    glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
    glyph3D->SetScaleModeToDataScalingOff();
    glyph3D->SetScaleFactor(surface->GetVtkPolyData()->GetLength()*scaleFactor);
    glyph3D->Update();

    //Create a mapper and actor for glyph
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(r,g,b);
    glyphActor->PickableOff();
    renderer->AddActor(glyphActor);
}

void FourChamberLandmarksView::PickCallBack() {

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(1E-4 * surface->GetVtkPolyData()->GetLength());
    int* eventPosition = interactor->GetEventPosition();
    int result = picker->Pick(float(eventPosition[0]), float(eventPosition[1]), 0.0, renderer);
    if (result == 0) return;
    double* pickPosition = picker->GetPickPosition();
    vtkIdList* pickedCellPointIds = surface->GetVtkPolyData()->GetCell(picker->GetCellId())->GetPointIds();

    double distance;
    int pickedSeedId = -1;
    double minDistance = 1E10;
    for (int i=0; i<pickedCellPointIds->GetNumberOfIds(); i++) {
        distance = vtkMath::Distance2BetweenPoints(
                    pickPosition, surface->GetVtkPolyData()->GetPoint(pickedCellPointIds->GetId(i)));
        if (distance < minDistance) {
            minDistance = distance;
            pickedSeedId = pickedCellPointIds->GetId(i);
        }//_if
    }//_for
    if (pickedSeedId == -1){
        pickedSeedId = pickedCellPointIds->GetId(0);
    }

    // update point
    pickedPoint.AddPointFromSurface(surface, pickedSeedId);

    m_Controls.widget_1->GetRenderWindow()->Render();
}

void FourChamberLandmarksView::KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*) {

    FourChamberLandmarksView* self;
    self = reinterpret_cast<FourChamberLandmarksView*>(ClientData);
    std::string key = self->interactor->GetKeySym();

    if (key == "space") {
        //Ask the labels
        self->PickCallBack();
        self->UserSelectLabel();

    } else if (key == "Delete") {

        self->pickedPoint.CleanupLastPoint();

        self->m_Controls.widget_1->GetRenderWindow()->Render();
    } else if (key == "H" || key == "h"){
        self->Help();
    }
}

// helper functions
void FourChamberLandmarksView::UserSelectLabel(){
    int dialogCode = inputsPickedPoints->exec();
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    int x = (screenGeometry.width() - inputsPickedPoints->width()) / 2;
    int y = (screenGeometry.height() - inputsPickedPoints->height()) / 2;
    inputsPickedPoints->move(x,y);

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        int pickedIndex = m_PickedPoints.comboBox->currentIndex();
        pickedPoint.PushBackLabelFromAvailable(pickedIndex);
        
    } else if (dialogCode == QDialog::Rejected) {
        inputsPickedPoints->close();
    }//_if
}

/*
========================
 CemrgApp radiobtn codes
========================
=== Landmarks ===
SVC_POST - radioBtn_RA_SVC_POST - 29
IVC_POST - radioBtn_RA_IVC_POST - 31
RAA_VALVE_P - radioBtn_RAA_TCV -  33
CS_VALVE_P - radioBtn_RA_CS_TCV - 35
SVC_ANT  - radioBtn_RA_SVC_ANT -  37
IVC_ANT  - radioBtn_RA_IVC_ANT -  39

=== Region ===
IVC_ANT - radioBtn_RA_IVC_ANT -      29
CS_TOP - radioBtn_RA_CS -            31
IVC_SVC_ANT - radioBtn_RA_IvcSvc -   33
SVC_ANT - radioBtn_RA_SVC_ANT -      35
RAA_ANT - radioBtn_RAA_ANT -         37
RAA_CS_ANT - radioBtn_RAA_CS_ANT -   39
*/
