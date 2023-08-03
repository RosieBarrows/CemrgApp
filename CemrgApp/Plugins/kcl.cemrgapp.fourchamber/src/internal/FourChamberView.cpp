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

//Micro services
#include <usModuleRegistry.h>
#ifdef _WIN32
// _WIN32 = we're in windows
#include <winsock2.h>
#else
// or linux/mac
#include <arpa/inet.h>
#endif

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>
#include <berryFileEditorInput.h>

//Qmitk
#include <mitkImage.h>
#include <mitkLog.h>
#include <QmitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkIDataStorageService.h>
// #include <mitkNodePredicateNot.h>
// #include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkLabelSetImage.h>

//VTK
#include <vtkFieldData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkTimerLog.h>
#include <vtkClipPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

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

// plugin classes
#include "kcl_cemrgapp_fourchamber_Activator.h"
#include "FourChamberView.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QDirIterator>
#include <QSignalMapper>
#include <QJsonObject>
#include <QJsonArray>
#include <QDate>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgFourChamberCmd.h>


const std::string FourChamberView::VIEW_ID = "org.mitk.views.fourchamberheart";

// JSON pre-determined filenames
const QString FourChamberView::POINTS_FILE = "physical_points.json";
const QString FourChamberView::POINTS_FILE_INDEX = "points.json";
const QString FourChamberView::GEOMETRY_FILE = "geometry.json";
const QStringList FourChamberView::SEGMENTATION_LIST = {"MIXED_LABEL", "BLOODPOOL", "LEFT_VENTRICLE", "RIGHT_VENTRICLE", "LEFT_ATRIUM", "RIGHT_ATRIUM", "AORTA", "LSPV", "LIPV", "RSPV", "RIPV", "LAA", "SVC", "IVC"};

void FourChamberView::SetFocus() {
    m_Controls.button_setfolder->setFocus();
}

void FourChamberView::CreateQtPartControl(QWidget *parent){
    // Link the slots to the buttons in the GUI here
    m_Controls.setupUi(parent);
    connect(m_Controls.button_setfolder, SIGNAL(clicked()), this, SLOT(SetWorkingFolder()));
    connect(m_Controls.button_prepseg, SIGNAL(clicked()), this, SLOT(PrepareSegmentation()));
    connect(m_Controls.button_meshing, SIGNAL(clicked()), this, SLOT(Meshing()));
    connect(m_Controls.button_uvcs, SIGNAL(clicked()), this, SLOT(UVCs()));
    connect(m_Controls.button_ventfibres, SIGNAL(clicked()), this, SLOT(VentricularFibres()));
    connect(m_Controls.button_atrfibres, SIGNAL(clicked()), this, SLOT(AtrialFibres()));
    connect(m_Controls.button_simset, SIGNAL(clicked()), this, SLOT(SimulationSetup()));

    connect(m_Controls.button_loaddicom, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_origin_spacing, SIGNAL(clicked()), this, SLOT(GetOriginSpacing()));
    connect(m_Controls.button_segment_image, SIGNAL(clicked()), this, SLOT(SegmentImgs()));
    connect(m_Controls.button_corrections, SIGNAL(clicked()), this, SLOT(Corrections()));
    connect(m_Controls.button_corrections_split, SIGNAL(clicked()), this, SLOT(CorrectionGetLabels()));
    connect(m_Controls.button_select_pts, SIGNAL(clicked()), this, SLOT(SelectPoints()));
    connect(m_Controls.button_select_pts_a, SIGNAL(clicked()), this, SLOT(SelectPointsCylinders()));
    connect(m_Controls.button_select_pts_b, SIGNAL(clicked()), this, SLOT(SelectPointsSlicers()));
    connect(m_Controls.button_select_pts_c, SIGNAL(clicked()), this, SLOT(SelectPointsValvePlains()));
    connect(m_Controls.button_select_pts_check, SIGNAL(clicked()), this, SLOT(SelectPointsCheck()));
    connect(m_Controls.button_select_pts_reset, SIGNAL(clicked()), this, SLOT(SelectPointsReset()));

    // Set default variables and initialise objects
    m_Controls.button_loaddicom->setVisible(false);
    m_Controls.button_origin_spacing->setVisible(false);
    m_Controls.button_segment_image->setVisible(false);

    m_Controls.button_select_pts->setVisible(false);
    m_Controls.button_select_pts_a->setVisible(false);
    m_Controls.button_select_pts_b->setVisible(false);
    m_Controls.button_select_pts_c->setVisible(false);
    m_Controls.button_select_pts_check->setVisible(false);
    m_Controls.button_select_pts_reset->setVisible(false);
    m_Controls.button_corrections->setVisible(false);
    m_Controls.button_corrections_split->setVisible(false);
    m_Controls.combo_corrections_id->setVisible(false);

    m_Controls.button_extractsurfs->setVisible(false);
    m_Controls.button_uvclandmarks->setVisible(false);
    m_Controls.button_calcuvcs->setVisible(false);

    SetButtonsEnableToOff();
    InitialiseJsonObjects();
    carpless = false;

    fourch_tools = std::unique_ptr<CemrgFourChamberTools>(new CemrgFourChamberTools());
    meshing_parameters = M3DParameters();
}

void FourChamberView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

// SLOTS
void FourChamberView::SetWorkingFolder(){
    if (!RequestProjectDirectoryFromUser()){
        MITK_WARN << "Folder not set. Check LOGFILE."; 
    } else{ // folder set successfully

        MITK_INFO << "Creating folder structure.";
        QStringList qsl = SDIR.Subdirectories();
        for (int ix=0; ix < qsl.size(); ix++) {
            MITK_INFO(QDir().mkpath(directory + "/" + qsl.at(ix))) << ("Folder created: [" + qsl.at(ix) + "]");
        }

        bool load_geometry_file = CheckForExistingFile(directory, FourChamberView::GEOMETRY_FILE);
        if (load_geometry_file) {
            MITK_INFO << "Loading geometry.json file";
            json_geometry = CemrgCommonUtils::ReadJSONFile(directory, FourChamberView::GEOMETRY_FILE);
            if (json_geometry["origin"].isUndefined() || json_geometry["spacing"].isUndefined()) {
                MITK_INFO << "Origin and/or spacing not defined in geometry.json file";
                m_Controls.button_origin_spacing->setEnabled(true);
            } else {
                ParseJsonArray(json_geometry, "origin", origin);
                ParseJsonArray(json_geometry, "spacing", spacing);

                std::cout << ArrayToString(origin, 3, "Loaded Origin").toStdString();
                std::cout << ArrayToString(spacing, 3, "Loaded Spacing").toStdString();
                m_Controls.button_origin_spacing->setEnabled(false);
            }
        }

        SetButtonsEnableToOn();
        m_Controls.button_loaddicom->setVisible(true);
        m_Controls.button_origin_spacing->setVisible(true);
        m_Controls.button_segment_image->setVisible(true);
    }
}

void FourChamberView::LoadDICOM() {
    int reply1 = QMessageBox::No;
#if defined(__APPLE__)
    MITK_INFO << "Ask user about alternative DICOM reader";
    reply1 = Ask("Question", "Use alternative DICOM reader?");
#endif

    if (reply1 == QMessageBox::Yes) {

        QString dicomFolder = QFileDialog::getExistingDirectory(NULL, "Open folder with DICOMs.", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        QString tmpNiftiFolder = cmd->DockerDicom2Nifti(dicomFolder);

        if (tmpNiftiFolder.compare("ERROR_IN_PROCESSING") != 0) {

            // add results in NIIs folder to Data Manager
            MITK_INFO << ("Conversion successful. Intermediate NII folder: " + tmpNiftiFolder).toStdString();
            Inform("Information", "Conversion successful, press the Process Images button to continue.");
            QDir niftiFolder(tmpNiftiFolder);
            QStringList niftiFiles = niftiFolder.entryList();

            if (niftiFiles.size()>0) {

                QString thisFile, path;
                for(int ix=0; ix<niftiFiles.size(); ix++) {

                    // load here files
                    thisFile = niftiFiles.at(ix);
                    if (thisFile.contains(".nii", Qt::CaseSensitive)) {
                        if (thisFile.contains("lge", Qt::CaseInsensitive) ||  thisFile.contains("mra", Qt::CaseInsensitive)) {

                            path = niftiFolder.absolutePath() + "/" + thisFile;
                            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
                            std::string key = "dicom.series.SeriesDescription";
                            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), thisFile.left(thisFile.length()-4).toStdString().c_str());

                        }//_if
                    }//_if
                }//_for

            } else {

                MITK_WARN << "Problem with conversion.";
                QMessageBox::warning(NULL, "Attention", "Problem with alternative conversion. Try MITK Dicom editor?");
                return;

            }//_if
        }//_if

    } else {

        MITK_INFO << "Using MITK DICOM editor";
        QString editor_id = "org.mitk.editors.dicomeditor";
        berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
        this->GetSite()->GetPage()->OpenEditor(input, editor_id);

    }//_if
}

void FourChamberView::GetOriginSpacing() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        Warn("Attention", "Please load and select only the image from the Data Manager to convert!");
        return;
    }//_if

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    mitk::BaseData::Pointer data = nodes[0]->GetData();
    if (data) {
        mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(data.GetPointer());
        if (image) {
            image->GetGeometry()->GetOrigin().ToArray(origin);
            image->GetGeometry()->GetSpacing().ToArray(spacing);

            QString origin_str = QString::number(origin[0]) + "," + QString::number(origin[1]) + "," + QString::number(origin[2]);
            QString spacing_str = QString::number(spacing[0]) + "," + QString::number(spacing[1]) + "," + QString::number(spacing[2]);

            QJsonObject geom_json = CemrgCommonUtils::CreateJSONObject({"origin", "spacing"}, {origin_str, spacing_str}, {"array", "array"});
            CemrgCommonUtils::WriteJSONFile(geom_json, directory, FourChamberView::GEOMETRY_FILE);

            MITK_INFO << ("Origin: (" + origin_str + ")").toStdString();
            MITK_INFO << ("Spacing: (" + spacing_str + ")").toStdString();
        }
    }
}

void FourChamberView::SegmentImgs() {
    int reply_load = Ask("Question", "Do you have a segmentation to load?");
    QString path = "";
    if (reply_load == QMessageBox::Yes) {
        path = QFileDialog::getOpenFileName(NULL, "Open Segmentation File", StdStringPath(SDIR.SEG).c_str(), QmitkIOUtil::GetFileOpenFilterString());

    } else if (reply_load == QMessageBox::No) {
        Inform("Attention", "Creating Multilabel Segmentation From CT data.\nSelect DICOM folder");
        QString dicom_folder = QFileDialog::getExistingDirectory(NULL, "Open DICOM folder",
             StdStringPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

        if (dicom_folder.isEmpty()) return;

        int reply_saveas = Ask("Question", "Do you want to save the DICOM as NIFTI?");
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);
        path = cmd->DockerCctaMultilabelSegmentation(Path(SDIR.SEG), dicom_folder, (reply_saveas==QMessageBox::Yes));
    }

    if (path.isEmpty()) return;

    double seg_origin[3], seg_spacing[3];
    mitk::Image::Pointer segmentation = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
        
    segmentation->GetVtkImageData()->GetSpacing(seg_spacing);
    segmentation->GetGeometry()->GetOrigin().ToArray(seg_origin);

    mitk::Point3D mitk_origin;
    mitk_origin[0] = origin[0];
    mitk_origin[1] = origin[1];
    mitk_origin[2] = origin[2];

    int test_origin = 0;
    for (int ix = 0; ix < 3; ix++){
        test_origin += (seg_origin[ix] == origin[ix]) ? 1:0;
    }
        
    if (!ArrayEqual(seg_origin, origin, 3)) {
        std::cout << test_origin << " Origin is different" << '\n';
        segmentation->SetOrigin(mitk_origin);
    }

    mitk::Vector3D mitk_spacing;
    mitk_spacing[0] = spacing[0];
    mitk_spacing[1] = spacing[1];
    mitk_spacing[2] = spacing[2];

    if (!ArrayEqual(seg_spacing, spacing, 3)){
        // spacing is different
        std::cout << "Spacing is different" << '\n';
        segmentation->SetSpacing(mitk_spacing);
    }

    QFileInfo fi(path);
    try {
        mitk::LabelSetImage::Pointer mlseg = mitk::LabelSetImage::New();
        mlseg->InitializeByLabeledImage(segmentation);
        mlseg->SetGeometry(segmentation->GetGeometry());

        CemrgCommonUtils::AddToStorage(mlseg, fi.baseName().toStdString(), this->GetDataStorage());
        mlseg->Modified();
    } catch (mitk::Exception &e) {
        MITK_ERROR << "Exception caught: " << e.GetDescription();
        CemrgCommonUtils::AddToStorage(segmentation, fi.baseName().toStdString(), this->GetDataStorage());
    }

    std::cout << ArrayToString(spacing, 3, "Spacing").toStdString();
    std::cout << ArrayToString(origin, 3, "Origin").toStdString() ;

    std::cout << ArrayToString(seg_spacing, 3, "SEG Spacing").toStdString();
    std::cout << ArrayToString(seg_origin, 3, "SEG Origin").toStdString();

}

void FourChamberView::PrepareSegmentation() {

    if (!RequestProjectDirectoryFromUser()) return;

    points_file_loaded = CheckForExistingFile(directory, FourChamberView::POINTS_FILE);
    if (points_file_loaded) {
        MITK_INFO << "Loading points.json file";
        json_points = CemrgCommonUtils::ReadJSONFile(directory, FourChamberView::POINTS_FILE);
        // iterate over json file keys. Unlock buttons if keys are all zeros
    } else {
        MITK_INFO << "Creating points.json file";
        CemrgCommonUtils::WriteJSONFile(json_points, directory, FourChamberView::POINTS_FILE);   
    }

    if (m_Controls.button_select_pts->isVisible()){

        m_Controls.button_corrections->setVisible(false);
        m_Controls.button_corrections_split->setVisible(false);
        m_Controls.button_select_pts->setVisible(false);
        m_Controls.button_select_pts_reset->setVisible(false);

    } else {
        m_Controls.button_corrections->setVisible(true);
        m_Controls.button_corrections_split->setVisible(true);
        m_Controls.button_select_pts->setVisible(true);
        m_Controls.button_select_pts_reset->setVisible(true);
    }

}

void FourChamberView::Meshing(){

    if (!RequestProjectDirectoryFromUser()) return;

    QString seg_dir = Path(SDIR.SEG);
    QString mesh_dir = Path(SDIR.MESH);
    QString segname = "seg_final_smooth.nrrd";

    if (!QDir().mkpath(mesh_dir)) {
        Warn("Problem with subfolder", "Problem with meshing subfolder"); 
        return;
    }

    int reply_load = Ask("Question", "Do you have a mesh to load?");
    QString mesh_path = "";

    // create cmd object (CemrgCommandLine) outputs to directory/meshing
    std::unique_ptr<CemrgCommandLine> cmd_object(new CemrgCommandLine());
    if (reply_load == QMessageBox::Yes) {
        mesh_path = QFileDialog::getOpenFileName(NULL, "Open Mesh Points File (.pts)", mesh_dir.toStdString().c_str(), tr("Points file (*.pts)"));
        QFileInfo fi(mesh_path);
        meshing_parameters.out_name = fi.baseName();
    } else if (reply_load == QMessageBox::No) {
    
        QString path = seg_dir + "/" + segname;
        bool ask_to_load = true;
        if(QFile::exists(path)) {
            std::string msg = "Load previously smoothed segmentation (" + segname.toStdString() + ")?";
            int reply = QMessageBox::question(NULL, "Question", msg.c_str());
            if(reply==QMessageBox::Yes){
                MITK_INFO << "Loading segmentation from file";
                ask_to_load = false;
            } 
        } 

        if (ask_to_load) { 
            // retrieve the path of the segmentation after user has selected it from the UI
            path = QFileDialog::getOpenFileName(NULL, "Open Segmentation file", seg_dir.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
            QFileInfo fi(path);
            seg_dir = fi.absolutePath();
            segname = fi.fileName();
        }

        QString inr_path = CemrgCommonUtils::ConvertToInr(seg_dir, segname, true, "converted.inr");
        if (inr_path.isEmpty()) {
            Warn("Failed to convert segmentation", "Error converting segmentation to INR format");
            return; 
        }

        QString parpath = mesh_dir + "/heart_mesh_data.par";
        bool success_m3d_params = UserSelectMeshtools3DParameters(inr_path);
        PrintMeshingParameters(parpath);

        if (success_m3d_params) {
            QString out_ext = meshing_parameters.out_carp ? "pts" : "vtk";

            QFileInfo finr(inr_path);
            Inform("Attention", "This operation takes a few minutes");
            QString mesh_path = cmd_object->ExecuteCreateCGALMesh(directory, meshing_parameters.out_name, parpath, SDIR.SEG + "/" + finr.fileName(), SDIR.MESH, out_ext);
        } 
    }

    MITK_INFO << ("Mesh path: " + mesh_path).toStdString();
    if (mesh_path.isEmpty()) return;

    // meshtool extract surface -msh=myocardium -ifmt=carp_txt -surf=test  -ofmt=vtk_polydata
    QStringList arguments;
    arguments << "extract"
              << "surface";
    arguments << "-msh=" + meshing_parameters.out_name;
    arguments << "-surf=" + meshing_parameters.out_name;
    arguments << "-ifmt=carp_txt";
    arguments << "-ofmt=vtk_polydata";

    QString outname = meshing_parameters.out_name + ".surfmesh";
    cmd_object->SetDockerImageOpenCarp();
    QString visualisation_surf_path = cmd_object->ExecuteCustomDocker(cmd_object->GetDockerImage(), mesh_dir, "meshtool", arguments, outname+".vtk");

    QFileInfo fvtk(visualisation_surf_path);

    visualisation_surf_path = cmd_object->DockerConvertMeshFormat(mesh_dir, outname, "vtk", outname, "vtk_polydata", 0.001);
    std::cout << "Visualisation surface loaded" << visualisation_surf_path << '\n';
    if (fvtk.exists()) {
        mitk::Surface::Pointer visualisation_surf = mitk::IOUtil::Load<mitk::Surface>(visualisation_surf_path.toStdString());

        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->Translate(origin[0], origin[1], origin[2]);

        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->SetTransform(transform);
        transformFilter->SetInputData(visualisation_surf->GetVtkPolyData());
        transformFilter->Update();

        visualisation_surf->SetVtkPolyData(transformFilter->GetOutput());

            CemrgCommonUtils::AddToStorage(visualisation_surf, fvtk.baseName().toStdString(), this->GetDataStorage());
    }
}

void FourChamberView::SelectLARALandmarks(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK!");
    }
}

void FourChamberView::ExtractSurfaces(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK!");
    }
}

void FourChamberView::SimulationSetup(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK");
    }
}

void FourChamberView::AtrialFibres(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK");
    }
}

void FourChamberView::VentricularFibres(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK");
    }
}

void FourChamberView::UVCs(){
    int reply = Ask("Question", "Is this the mesh for which you would like to calculate UVCs?");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK!");
    }
    if(reply==QMessageBox::No){
        Inform("Answer", "Choose a different mesh.");
    }

    if (m_Controls.button_extractsurfs->isVisible()) {
        m_Controls.button_extractsurfs->setVisible(true);
    } else {
        m_Controls.button_extractsurfs->setVisible(true);
    }

    if (m_Controls.button_uvclandmarks->isVisible()) {
        m_Controls.button_uvclandmarks->setVisible(true);
    } else {
        m_Controls.button_uvclandmarks->setVisible(true);
    }

    if (m_Controls.button_calcuvcs->isVisible()) {
        m_Controls.button_calcuvcs->setVisible(true);
    } else {
        m_Controls.button_calcuvcs->setVisible(true);
    }
}

void FourChamberView::CalculateUVCs(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        Inform("Answer", "OK!");
    }
}

void FourChamberView::Corrections(){

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    bool button_timing = m_Controls.button_corrections->text().contains("PV Split");
    
    std::string msg = "Loading Multilabel Segmentation tools.\n\n";
    msg += button_timing ? " Make any manual corrections on pulmonary veins." : "Make additional corrections to segmentation";
    
    Inform("Attention", msg.c_str());
    
    if (button_timing) {
        m_Controls.button_corrections->setText("        2.1: Manual Corrections");
        m_Controls.button_corrections->setStyleSheet("QPushButton {text-align: left; color: rgb(132, 174, 235);}");
    }

    this->GetSite()->GetPage()->ShowView("org.mitk.views.multilabelsegmentation");
}

void FourChamberView::CorrectionGetLabels() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        Warn("Attention", "Please load and select only the Segmentation from the Data Manager to convert!");
        return;
    }//_if

    mitk::BaseData::Pointer data = nodes[0]->GetData();
    std::string nodeName = nodes[0]->GetName();
    if (data) {
        mitk::Image::Pointer seg = dynamic_cast<mitk::Image *>(data.GetPointer());
        if (seg) {
            
            // find which labels need splitting
            fourch_tools->GetLabels(seg, labelsInSegmentation);
           
            mitk::LabelSetImage::Pointer mlseg = mitk::LabelSetImage::New();
            mlseg->InitializeByLabeledImage(seg);
            mlseg->SetGeometry(seg->GetGeometry());

            foreach (int label, labelsInSegmentation) {
                QString cbox_label = QString::number(label);
                m_Controls.combo_corrections_id->addItem(QString::number(label));

                mlseg->GetLabel(label)->SetOpacity(0.3);
            }
            mlseg->Modified();
            CemrgCommonUtils::UpdateFromStorage(mlseg, nodeName, this->GetDataStorage());

            m_Controls.combo_corrections_id->setVisible(true);
            connect(m_Controls.combo_corrections_id, SIGNAL(activated(int)), this, SLOT(CorrectionIdLabels(int)));
        } // _if_image
    } // _if_data
}

void FourChamberView::CorrectionIdLabels(int index) {
    if (m_Controls.combo_corrections_id->count() == 0)
        return;
    
     //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        Warn("Attention", "Please load and select only the Segmentation from the Data Manager to convert!");
        return;
    }//_if

    QString cbox_label = m_Controls.combo_corrections_id->itemText(index);
    unsigned int label = cbox_label.split(":")[0].toInt();
    mitk::BaseData::Pointer data = nodes[0]->GetData();
    if (data) {
        mitk::LabelSetImage::Pointer mlseg = dynamic_cast<mitk::LabelSetImage *>(data.GetPointer());
        if (mlseg) {
            mlseg->GetLabel(label)->SetOpacity(1.0);
            mlseg->Modified();
            CemrgCommonUtils::UpdateFromStorage(mlseg, nodes[0]->GetName(), this->GetDataStorage());

            QColor qcolor = MitkColorToQColor(mlseg->GetLabel(label)->GetColor());
            bool userInputsAccepted = UserSelectIdentifyLabels(label, qcolor);
            if (userInputsAccepted) {
                if (splitCurrentLabel) {
                    mitk::Image::Pointer updateSeg = fourch_tools->SplitLabelsOnRepeat((mitk::Image::Pointer) mlseg, (int) label);
                    mlseg->InitializeByLabeledImage(updateSeg);
                    mlseg->SetGeometry(updateSeg->GetGeometry());

                    CemrgCommonUtils::UpdateFromStorage(mlseg, nodes[0]->GetName(), this->GetDataStorage());
                } else {
                    m_Controls.combo_corrections_id->setItemText(index, QString::number(label) + ": Identified");
                }
            }
            MITK_INFO(userInputsAccepted) << "User inputs accepted";
        } // _if_image
    }
}

void FourChamberView::SelectPoints() {

    if (!RequestProjectDirectoryFromUser()) return;

    if (!m_Controls.button_select_pts_a->isVisible()) {
        m_Controls.button_select_pts_a->setVisible(true);
        m_Controls.button_select_pts_b->setVisible(true);
        m_Controls.button_select_pts_c->setVisible(true);
        m_Controls.button_select_pts_check->setVisible(true);
    } else {
        m_Controls.button_select_pts_a->setVisible(false);
        m_Controls.button_select_pts_b->setVisible(false);
        m_Controls.button_select_pts_c->setVisible(false);
        m_Controls.button_select_pts_check->setVisible(false);
        return; 
    }

    std::string pset_name = "four_chamber_point_set";

    // Create pointset
    pset = mitk::PointSet::New();
    CemrgCommonUtils::AddToStorage(pset, pset_name, this->GetDataStorage(), false);

    Inform("Controls", "SHIFT + Left Click to add points");

    CemrgDataInteractor::Pointer m_interactor = CemrgDataInteractor::New();
    m_interactor->Initialise(directory + "/" + FourChamberView::POINTS_FILE);
    m_interactor->LoadStateMachine("PointSet.xml");
    m_interactor->SetEventConfig("PointSetConfig.xml");
    m_interactor->SetDataNode(this->GetDataStorage()->GetNamedNode(pset_name.c_str()));

    if (points_file_loaded) { 
        // load points from json file
        QStringList allKeys = json_points.keys();
        foreach (QString Key, allKeys) {
            double arr[3] = {0,0,0};
            ParseJsonArray(json_points, Key, arr);
            pset->InsertPoint(atoi(Key.toStdString().c_str()), mitk::Point3D(arr));
        }
    }
}

void FourChamberView::SelectPointsCylinders() {
    
    if (m_Controls.button_select_pts_a->text().contains("Check Cylinder")) {
        if (!CheckPointsInJsonFile(ManualPointsType::CYLINDERS)) {
            Warn("Attention - Missing Points", "All points need to be defined before running the script.");
            SelectPointsCheck();
            return;
        }
        m_Controls.button_select_pts_a->setText("Run Scripts for Cylinders");
    } else {
        // Run cylinders script
        Inform("Attention", "Running Scripts for Cylinders");
        m_Controls.button_select_pts_a->setEnabled(false);

    }
}

void FourChamberView::SelectPointsSlicers() {
    
    if (m_Controls.button_select_pts_b->text().contains("Check Slicer")) {
        if (!CheckPointsInJsonFile(ManualPointsType::SLICERS)) {
            Warn("Attention - Missing Points", "All points need to be defined before running the script.");
            SelectPointsCheck();
            return;
        }
        m_Controls.button_select_pts_b->setText("Run Scripts for Slicers");
    } else {
        // Run Slicers script
        Inform("Attention", "Running Scripts for Slicers");
        m_Controls.button_select_pts_b->setEnabled(false);
    }
}

void FourChamberView::SelectPointsValvePlains(){
    
    if (m_Controls.button_select_pts_c->text().contains("Check Valve Plains")) {
        if (!CheckPointsInJsonFile(ManualPointsType::VALVE_PLAINS)) {
            Warn("Attention - Missing Points", "All points need to be defined before running the script.");
            SelectPointsCheck();
            return;
        }
        m_Controls.button_select_pts_c->setText("Run Scripts for Valve Plains");
    } else {
        // Run Valve Plains script
        Inform("Attention", "Running Scripts for Valve Plains");
        m_Controls.button_select_pts_c->setEnabled(false);
    }
}

bool FourChamberView::CheckPointsInJsonFile(ManualPointsType mpt){
    // reload json file
    ReloadJsonPoints();

    QStringList pointsToCheck = SegPointIds.GetPointLabelOptions(mpt);
    bool allPointsExist = true;
    double testArray[3] = {0,0,0};
    foreach (QString Key, pointsToCheck) {
        double arr[3] = {0,0,0};
        ParseJsonArray(json_points, Key, arr);
        if (ArrayEqual(arr, testArray, 3)) {
            std::cout << ArrayToString(arr, 3, Key).toStdString();
            return false;
        }
    }
    return allPointsExist;
}

void FourChamberView::ReloadJsonPoints(){
    QJsonObject json = CemrgCommonUtils::ReadJSONFile(directory, FourChamberView::POINTS_FILE);
    QStringList allKeys = json_points.keys();

    int countDifferent = 0;
    foreach (QString Key, allKeys) {
        double arrCurrent[3] = {0,0,0};
        double arrNew[3] = {0,0,0};
        ParseJsonArray(json, Key, arrNew);
        ParseJsonArray(json_points, Key, arrCurrent);
        if (!ArrayEqual(arrCurrent, arrNew, 3)) {
            countDifferent++;
        }
    }

    if (countDifferent > 0) {
        // This assumes that the file is always newer than what the points are already saved.
        json_points = json;
    }
}

void FourChamberView::SelectPointsReset(){
    MITK_INFO << "Resetting points";
    int reply = Ask("Question", "Are you sure you want to reset the points?\nThis will overwrite the file too.");
    if(reply==QMessageBox::Yes){
        InitialiseJsonObjects();
        CemrgCommonUtils::WriteJSONFile(json_points, directory, FourChamberView::POINTS_FILE);

        // remove pointset from storage
        this->GetDataStorage()->Remove(this->GetDataStorage()->GetNamedNode("four_chamber_point_set"));

        // reset buttons
        m_Controls.button_select_pts_a->setText("Check Cylinders Points");
        m_Controls.button_select_pts_a->setEnabled(true);
        m_Controls.button_select_pts_b->setText("Check Slicers Points");
        m_Controls.button_select_pts_b->setEnabled(true);
        m_Controls.button_select_pts_c->setText("Check Valve Plains Points");
        m_Controls.button_select_pts_c->setEnabled(true);
    }

}

void FourChamberView::SelectPointsCheck() {
    
    QJsonObject points = CemrgCommonUtils::ReadJSONFile(directory, FourChamberView::POINTS_FILE);
    QStringList cylinders = SegPointIds.CYLINDERS();
    QStringList slicers = SegPointIds.SLICERS();
    QStringList valve_plains = SegPointIds.VALVE_PLAINS();

    std::string output = "Selected Points\n\n";

    output += PrintPoints(points, cylinders, "Cylinders") + '\n';
    output += PrintPoints(points, slicers, "Slicers") + '\n';
    output += PrintPoints(points, valve_plains, "Valve Plains") + '\n';

    std::cout << output.c_str();
    QMessageBox::information(nullptr, "Selected Points", output.c_str());
}

std::string FourChamberView::PrintPoints(QJsonObject json, QStringList keysList, QString title) {
    std::string output = title.toStdString() + "\n";
    for (int ix = 0; ix < keysList.size(); ix++) {
        double arr[3] = {0, 0, 0};
        QString key = keysList.at(ix);
        if (json[key].isUndefined()) {
            output += key.toStdString() + ": Undefined\n";
            continue;
        }
        output += key.toStdString() + ": ";
        for (int jx = 0; jx < 3; jx++) {
            arr[jx] = json[key].toArray().at(jx).toDouble();
        }

        if (arr[0] == 0 && arr[1] == 0 && arr[2] == 0) {
            output += "NOT SET\n";
            continue;
        }

        for (int jx = 0; jx < 3; jx++) {
            output += QString::number(arr[jx]).toStdString();
            output += (jx < 2) ? ", " : ")\n";
        }
    }

    return output;
}

// helper`
bool FourChamberView::RequestAnyFolderFromUser(QString & dir, std::string msg, bool project_dir){

    bool successfulAssignment = true;

    // ask the user for a directory 
    if (dir.isEmpty()) {
        
        MITK_INFO << "Directory is empty. Requesting user for dirrectory.";
        dir = QFileDialog::getExistingDirectory( NULL, msg.c_str(),
            mitk::IOUtil::GetProgramPath().c_str(),QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

        MITK_INFO << ("Directory selected:" + dir).toStdString();

        if (dir.isEmpty() || dir.simplified().contains(" ")) {
            MITK_WARN << "Please select a project directory with no spaces in the path!";
            QMessageBox::warning(NULL, "Attention", "Please select a project dir with no spaces in the path!");
            dir = QString();
            successfulAssignment = false;
        }//_if

        if (successfulAssignment && project_dir){
            QString now = QDate::currentDate().toString(Qt::ISODate);
            QString logfilename = dir + "/LOGFILE_" + now + ".log";
            std::string logfname_str = logfilename.toStdString();

            mitk::LoggingBackend::SetLogFile(logfname_str.c_str());
            MITK_INFO << ("Changed logfile location to: " + logfilename).toStdString();
        }

    } else {
        MITK_INFO << ("Folder already set: " + dir).toStdString();
    }//_if

    return successfulAssignment;

}

bool FourChamberView::CheckForExistingFile(QString dir, QString filename) {
    MITK_INFO << ("Checking for existing " + filename + " file").toStdString();
    QFileInfo fi(dir + "/" + filename);
    bool load_file = false;
    if (fi.exists()){
        std::string msg2 = "File " + filename.toStdString() + " already exists in: \n[" + dir.toStdString() + "]\n\nYES=Load it \t NO=Overwrite it";
        int reply_load_json = Ask( "Attention" , msg2.c_str());
        if (reply_load_json == QMessageBox::Yes){
            load_file = true;
        } else {
            MITK_INFO << ("Overwriting " + filename + " file").toStdString();
        }
    }
    return load_file;
}

int FourChamberView::Ask(std::string title, std::string msg){
    MITK_INFO << "[USER_INPUT]" + msg;
    return QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);
}

void FourChamberView::Warn(std::string title, std::string msg){
    MITK_WARN << msg;
    QMessageBox::warning(NULL, title.c_str(), msg.c_str());
}

void FourChamberView::Inform(std::string title, std::string msg){
    MITK_INFO << msg;
    QMessageBox::information(NULL, title.c_str(), msg.c_str());
}

void FourChamberView::InitialiseJsonObjects() {

    QStringList pt_keys = SegPointIds.GetPointLabelOptions(ManualPointsType::CYLINDERS);
    pt_keys.append(SegPointIds.GetPointLabelOptions(ManualPointsType::SLICERS));
    pt_keys.append(SegPointIds.GetPointLabelOptions(ManualPointsType::VALVE_PLAINS));

    QStringList values = QStringList(), types = QStringList();
    InitialiseQStringListsFromSize(pt_keys.size(), values, types);

    json_points = CemrgCommonUtils::CreateJSONObject(pt_keys, values, types);
    
}

void FourChamberView::InitialiseQStringListsFromSize(int num, QStringList& values, QStringList& types) {
    for (int ix = 0; ix < num; ix++) {
        values << "0.0,0.0,0.0";
        types << "array";
    }
}

void FourChamberView::PrintMeshingParameters(QString path_to_par){
    std::ofstream fo(path_to_par.toStdString());
    fo << "#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -" << '\n';
    fo << "#Data file for meshtools3d utility"<< '\n';
    fo << "#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -"<< '\n';
    fo << "[segmentation]"<< '\n';
    fo << "seg_dir = " << meshing_parameters.seg_dir.toStdString() << '\n';
    fo << "seg_name = " << meshing_parameters.seg_name.toStdString() << '\n';
    fo << "mesh_from_segmentation = " << meshing_parameters.mesh_from_segmentation << '\n';
    fo << "boundary_relabelling = " << meshing_parameters.boundary_relabelling << '\n';
    fo << "[meshing]"<< '\n';
    fo << "facet_angle         = "<< meshing_parameters.facet_angle << '\n';
    fo << "facet_size          = "<< meshing_parameters.facet_size << '\n';
    fo << "facet_distance      = "<< meshing_parameters.facet_distance << '\n';
    fo << "cell_rad_edge_ratio = "<< meshing_parameters.cell_rad_edge_ratio << '\n';
    fo << "cell_size           = "<< meshing_parameters.cell_size << '\n';
    fo << "rescaleFactor       = "<< meshing_parameters.rescaleFactor << '\n';
    fo << "[laplacesolver]" << '\n';
    fo << "abs_toll        = " << meshing_parameters.abs_tol << '\n';
    fo << "rel_toll        = " << meshing_parameters.rel_tol << '\n';
    fo << "itr_max         = " << meshing_parameters.itr_max << '\n';
    fo << "dimKrilovSp     = " << meshing_parameters.dimKrilovSp << '\n';
    fo << "verbose         = " << meshing_parameters.verbose << '\n';
    fo << "[others]"<< '\n';
    fo << "eval_thickness  = " << meshing_parameters.eval_thickness << '\n';
    fo << "[output]"<< '\n';
    fo << "outdir          = " << meshing_parameters.out_dir.toStdString() << '\n';
    fo << "name            = " << meshing_parameters.out_name.toStdString() << '\n';
    fo << "out_medit       = " << meshing_parameters.out_medit << '\n';
    fo << "out_carp        = " << meshing_parameters.out_carp << '\n';
    fo << "out_carp_binary = " << meshing_parameters.out_carp_binary << '\n';
    fo << "out_vtk         = " << meshing_parameters.out_vtk << '\n';
    fo << "out_vtk_binary  = " << meshing_parameters.out_vtk_binary << '\n';
    fo << "out_potential   = " << meshing_parameters.out_potential << '\n';
    fo.close();

    QStringList keys = QStringList();
    QStringList values = QStringList();
    QStringList types = QStringList();

    QFileInfo fi(path_to_par);
    meshing_parameters.KeysAndValues(keys, values, types);
    QJsonObject json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
    MITK_INFO(CemrgCommonUtils::WriteJSONFile(json, fi.absolutePath(), fi.baseName() + ".json"));
}

void FourChamberView::LoadMeshingParametersFromJson(QString dir, QString json_file){
    QJsonObject json = CemrgCommonUtils::ReadJSONFile(dir, json_file);
    meshing_parameters.seg_dir = json["seg_dir"].toString();
    meshing_parameters.seg_name = json["seg_name"].toString();
    meshing_parameters.mesh_from_segmentation = json["mesh_from_segmentation"].toBool();
    meshing_parameters.boundary_relabelling = json["boundary_relabelling"].toBool();
    meshing_parameters.facet_angle = json["facet_angle"].toDouble();
    meshing_parameters.facet_size = json["facet_size"].toDouble();
    meshing_parameters.facet_distance = json["facet_distance"].toDouble();
    meshing_parameters.cell_rad_edge_ratio = json["cell_rad_edge_ratio"].toDouble();
    meshing_parameters.cell_size = json["cell_size"].toDouble();
    meshing_parameters.rescaleFactor = json["rescaleFactor"].toDouble();
    meshing_parameters.abs_tol = json["abs_tol"].toDouble();
    meshing_parameters.rel_tol = json["rel_tol"].toDouble();
    meshing_parameters.itr_max = json["itr_max"].toInt();
    meshing_parameters.dimKrilovSp = json["dimKrilovSp"].toInt();
    meshing_parameters.verbose = json["verbose"].toBool();
    meshing_parameters.eval_thickness = json["eval_thickness"].toBool();
    meshing_parameters.out_dir = json["out_dir"].toString();
    meshing_parameters.out_name = json["out_name"].toString();
    meshing_parameters.out_medit = json["out_medit"].toBool();
    meshing_parameters.out_carp = json["out_carp"].toBool();
    meshing_parameters.out_carp_binary = json["out_carp_binary"].toBool();
    meshing_parameters.out_vtk = json["out_vtk"].toBool();
    meshing_parameters.out_vtk_binary = json["out_vtk_binary"].toBool();
    meshing_parameters.out_potential = json["out_potential"].toBool();
}

QString FourChamberView::GetPointTypeString(ManualPointsType mpt) {
    QString res;
    switch (mpt) {
        case ManualPointsType::CYLINDERS :
            res = "Cylinders";
            break;
        case ManualPointsType::SLICERS :
            res = "Slicers";
            break;
        case ManualPointsType::VALVE_PLAINS :
            res = "ValvePlains";
            break;
        default:
            res = "";
            break;
    }

    return res;
}

// User Select Functions 
bool FourChamberView::UserSelectMeshtools3DParameters(QString pre_input_path){
    bool userAccepted = false;
    QDialog* inputs = new QDialog(0, 0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_m3d.setupUi(inputs);
    connect(m_m3d.dialog_button, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_m3d.dialog_button, SIGNAL(rejected()), inputs, SLOT(reject()));

    connect(m_m3d.button_load_image, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_m3d.button_load_image, directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(M3dBrowseFile(const QString&)));

    QFileInfo fi(pre_input_path);
    m_m3d.lineEdit_input_path->setPlaceholderText(pre_input_path);
    m_m3d.lineEdit_out_dir->setText(directory + "/" + SDIR.MESH);

    int dialogCode = inputs->exec();

    // Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        // parse all booleans
        meshing_parameters.eval_thickness = m_m3d.checkBox_eval_thickness->isChecked();

        meshing_parameters.out_medit = m_m3d.checkBox_out_medit->isChecked();
        meshing_parameters.out_carp = m_m3d.checkBox_out_carp->isChecked();
        meshing_parameters.out_vtk = m_m3d.checkBox_out_vtk->isChecked();
    
        if (meshing_parameters.CheckFormats()){
            Inform("Attention", "No output format selected.");
            return false;
        }

        if (m_m3d.checkBox_out_binary->isChecked()) {
            meshing_parameters.SetBinaryFormats();
        }

        // parse all strings
        if (!m_m3d.lineEdit_input_path->text().isEmpty()) {
            fi.setFile(m_m3d.lineEdit_input_path->text());
        }
        meshing_parameters.seg_dir = fi.absolutePath();
        meshing_parameters.seg_name = fi.fileName();

        meshing_parameters.out_dir = m_m3d.lineEdit_out_dir->text();
        if (m_m3d.lineEdit_out_name->text().isEmpty()) {
            meshing_parameters.out_name = "myocardium";
        } else {
            meshing_parameters.out_name = m_m3d.lineEdit_out_name->text();
        }

        // parse all numbers
        bool ok1, ok2, ok3, ok4, ok5, ok6, ok7, ok8, ok9, ok10;
        meshing_parameters.facet_angle = m_m3d.lineEdit_facet_angle->text().toDouble(&ok1);
        meshing_parameters.facet_size = m_m3d.lineEdit_facet_size->text().toDouble(&ok3);
        meshing_parameters.facet_distance = m_m3d.lineEdit_facet_dist->text().toDouble(&ok2);
        meshing_parameters.cell_rad_edge_ratio = m_m3d.lineEdit_cell_rad_edge_ratio->text().toDouble(&ok4);

        meshing_parameters.cell_size = m_m3d.lineEdit_cell_size->text().toDouble(&ok5);
        meshing_parameters.rescaleFactor = m_m3d.lineEdit_rescaleFactor->text().toDouble(&ok6);

        meshing_parameters.itr_max = m_m3d.lineEdit_itr_max->text().toInt(&ok9);
        meshing_parameters.abs_tol = m_m3d.lineEdit_abs_tol->text().toDouble(&ok7);
        meshing_parameters.rel_tol = m_m3d.lineEdit_rel_tol->text().toDouble(&ok8);
        meshing_parameters.dimKrilovSp = m_m3d.lineEdit_dimKrilovSp->text().toInt(&ok10);

        if (!ok1) meshing_parameters.facet_angle = 30.0;
        if (!ok3) meshing_parameters.facet_size = 0.8;
        if (!ok2) meshing_parameters.facet_distance = 4;
        if (!ok4) meshing_parameters.cell_rad_edge_ratio = 2.0;
        if (!ok5) meshing_parameters.cell_size = 0.8;
        if (!ok6) meshing_parameters.rescaleFactor = 1000;
        if (!ok9) meshing_parameters.itr_max = 700;
        if (!ok7) meshing_parameters.abs_tol = 1e-6;
        if (!ok8) meshing_parameters.rel_tol = 1e-6;
        if (!ok10) meshing_parameters.dimKrilovSp = 500;

        userAccepted = true;
    }

    return userAccepted;
}

bool FourChamberView::UserSelectIdentifyLabels(unsigned int label, QColor qc){
    QDialog *inputs = new QDialog(0, 0);
    bool userInputAccepted = false;
    m_IdLabels.setupUi(inputs);

    m_IdLabels.label_colour->setStyleSheet("background-color: " + qc.name() + ";");
    m_IdLabels.combo_id_label->addItems(FourChamberView::SEGMENTATION_LIST);
    m_IdLabels.lineEdit_current_label->setText(QString::number(label));
    
    connect(m_IdLabels.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_IdLabels.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    
    int dialogCode = inputs->exec();
    if (dialogCode == QDialog::Accepted) {
        bool ok;
        imageLabel = label;
        defaultLabel = m_IdLabels.combo_id_label->currentIndex();
        userLabel = m_IdLabels.lineEdit_current_label->text().toInt(&ok);
        if (!ok) {
            userLabel = defaultLabel;
        }

        splitCurrentLabel = (m_IdLabels.combo_id_label->currentText() == "MIXED_LABEL");
        if (splitCurrentLabel) {
            Inform("Attention: Mixed Label", "Will attempt to split");
            defaultLabel = -1;
            userLabel = -1;
        }

        userInputAccepted = true;
    } else {
        inputs->deleteLater();
    }
    return userInputAccepted;
}

void FourChamberView::M3dBrowseFile(const QString &dir) {
    QString titlelabel, input = "";
    std::string msg;

    msg = "Select segmentation image for meshing";

    Inform("Attention", msg.c_str());
    input = QFileDialog::getOpenFileName(NULL, msg.c_str(), dir, QmitkIOUtil::GetFileOpenFilterString());

    QFileInfo fi(input);
    m_m3d.lineEdit_input_path->setText(input);
}

void FourChamberView::SetButtonsEnable(bool enable) {
    // Sets all buttons to enable or disable.
    // Used at the beginning of the pipeline
    m_Controls.button_prepseg->setEnabled(enable) ; 
    m_Controls.button_meshing->setEnabled(enable) ; 
    m_Controls.button_uvcs->setEnabled(enable) ; 
    m_Controls.button_ventfibres->setEnabled(enable) ; 
    m_Controls.button_atrfibres->setEnabled(enable) ; 
    m_Controls.button_simset->setEnabled(enable) ;
}

QString FourChamberView::ArrayToString(double *arr, int size, const QString &title) {
    QString msg =  title + ": (";
    for (int ix = 0; ix < size; ix++) {
        QString endstr = (ix < size - 1) ? "," : ")";
        msg += QString::number(arr[ix]); 
        msg += endstr;
    }
    msg += '\n';
    return msg;
}

QString FourChamberView::MitkColorToHex(const mitk::Color& colour) {
    QColor qcolour = MitkColorToQColor(colour);
        
    return qcolour.name();
}

QColor FourChamberView::MitkColorToQColor(const mitk::Color& colour) {
    int red = static_cast<int>(colour[0] * 255);
    int green = static_cast<int>(colour[1] * 255);
    int blue = static_cast<int>(colour[2] * 255);

    QColor res(red, green, blue);
    return res;
}

bool FourChamberView::ArrayEqual(double *arr1, double *arr2, int size, double tol) {
    for (int ix = 0; ix < size; ix++) {
        if (fabs(arr1[ix] - arr2[ix]) > tol) {
            return false;
        }
    }
    return true;
}

void FourChamberView::ParseJsonArray(QJsonObject json, QString key, double *arr, int size) {
    if (json[key].isUndefined()){
        MITK_WARN << ("Key [" + key + "] is Undefined for JSON object").toStdString();
    }
    for (int ix = 0; ix < size; ix++) {
        arr[ix] = json[key].toArray().at(ix).toDouble();
    }
}