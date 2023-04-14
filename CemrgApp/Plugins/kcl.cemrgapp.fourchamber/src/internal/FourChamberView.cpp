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
#include "CemrgDataInteractor.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QDirIterator>
#include <QSignalMapper>
#include <QDate>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>


const std::string FourChamberView::VIEW_ID = "org.mitk.views.fourchamberheart";
// JSON pre-determined filenames
const QString FourChamberView::POINTS_FILE = "points.json";
const QString FourChamberView::GEOMETRY_FILE = "geometry.json";

void FourChamberView::SetFocus(){
    // m_Controls.button_setfolder->setFocus();
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
    connect(m_Controls.button_splitpulms, SIGNAL(clicked()), this, SLOT(SplitPulmVeins()));
    connect(m_Controls.button_select_pts_a, SIGNAL(clicked()), this, SLOT(SelectPointsCylinders()));
    connect(m_Controls.button_select_pts_b, SIGNAL(clicked()), this, SLOT(SelectPointsSlicers()));
    connect(m_Controls.button_select_pts_c, SIGNAL(clicked()), this, SLOT(SelectPointsValvePlains()));
    connect(m_Controls.button_corrections, SIGNAL(clicked()), this, SLOT(Corrections()));

    // Set default variables and initialise objects
    m_Controls.button_loaddicom->setVisible(false);
    m_Controls.button_origin_spacing->setVisible(false);
    m_Controls.button_segment_image->setVisible(false);

    m_Controls.button_splitpulms->setVisible(false);
    m_Controls.button_select_pts_a->setVisible(false);
    m_Controls.button_select_pts_b->setVisible(false);
    m_Controls.button_select_pts_c->setVisible(false);
    m_Controls.button_corrections->setVisible(false);

    m_Controls.button_extractsurfs->setVisible(false);
    m_Controls.button_uvclandmarks->setVisible(false);
    m_Controls.button_calcuvcs->setVisible(false);


    InitialiseJsonObjects();
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

        MITK_INFO << "Checking for existing points.json file";
        QFileInfo fi(directory + "/" + FourChamberView::POINTS_FILE);

        if (fi.exists()){
            int reply_load_json = Ask( "Attention" ,"Load previously found points file?");
            if (reply_load_json == QMessageBox::Yes){
                json_points = CemrgCommonUtils::ReadJSONFile(directory, FourChamberView::POINTS_FILE);
            } else {
                QMessageBox::information(NULL, "Attention", "Overwriting points.json file");
            }
        }

        points_file_loaded = CemrgCommonUtils::WriteJSONFile(json_points, directory, FourChamberView::POINTS_FILE);
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
            MITK_INFO << ("Conversion succesful. Intermediate NII folder: " + tmpNiftiFolder).toStdString();
            QMessageBox::information(NULL, "Information", "Conversion successful, press the Process Images button to continue.");
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
        Warn("Attention", "Please load and select only the DICOM image from the Data Manager to convert!");
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
        QMessageBox::information(NULL, "Attention", "Creating Multilabel Segmentation From CT data.\nSelect DICOM folder");
        QString dicom_folder = QFileDialog::getExistingDirectory(NULL, "Open DICOM folder",
             StdStringPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

        if (dicom_folder.isEmpty()) return;

        int reply_saveas = Ask("Question", "Do you want to save the DICOM as NIFTI?");
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);
        path = cmd->DockerCctaMultilabelSegmentation(Path(SDIR.SEG), dicom_folder, (reply_saveas==QMessageBox::Yes));
    }

    if (path.isEmpty()) return;

    QFileInfo fi(path);
    double seg_origin[3], seg_spacing[3];
    mitk::Image::Pointer im = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
    im->GetVtkImageData()->GetSpacing(seg_spacing);
    im->GetGeometry()->GetOrigin().ToArray(seg_origin);

    std::cout << "Spacing: [" << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << "]\n";
    std::cout << "Origin: [" << origin[0] << ", " << origin[1] << ", " << origin[2] << "]\n";

    std::cout << "SEG Spacing: [" << seg_spacing[0] << ", " << seg_spacing[1] << ", " << seg_spacing[2] << "]\n";
    std::cout << "SEG Origin: [" << seg_origin[0] << ", " << seg_origin[1] << ", " << seg_origin[2] << "]\n";

    im->GetGeometry()->SetSpacing(spacing);
    im->GetGeometry()->SetOrigin(origin);

    try {
        mitk::LabelSetImage::Pointer mlseg = mitk::LabelSetImage::New();
        mlseg->InitializeByLabeledImage(im);
        CemrgCommonUtils::AddToStorage(mlseg, fi.baseName().toStdString(), this->GetDataStorage());
        mlseg->Modified();
    } catch (mitk::Exception &e) {
        MITK_ERROR << "Exception caught: " << e.GetDescription();
        CemrgCommonUtils::AddToStorage(im, fi.baseName().toStdString(), this->GetDataStorage());
    }

}

void FourChamberView::PrepareSegmentation() {

    if (!RequestProjectDirectoryFromUser()) return;

    if (m_Controls.button_select_pts_a->isVisible()){

        m_Controls.button_select_pts_a->setVisible(false);
        m_Controls.button_select_pts_b->setVisible(false);
        m_Controls.button_select_pts_c->setVisible(false);

    } else {
        m_Controls.button_splitpulms->setVisible(true);
        m_Controls.button_select_pts_a->setVisible(true);
        m_Controls.button_select_pts_b->setVisible(true);
        m_Controls.button_select_pts_c->setVisible(true);
        m_Controls.button_corrections->setVisible(true);
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
        mesh_path = QFileDialog::getOpenFileName(NULL, "Open Mesh Points File (.pts)", mesh_dir.toStdString().c_str(), tr("Parameter file (*.pts)"));
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
            QMessageBox::information(NULL, "Attention", "This operation takes a few minutes");
            QString mesh_path = cmd_object->ExecuteCreateCGALMesh(directory, meshing_parameters.out_name, parpath, SDIR.SEG + "/" + finr.fileName(), SDIR.MESH, out_ext);
        } 
    }

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
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::ExtractSurfaces(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::SimulationSetup(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK");
    }
}

void FourChamberView::AtrialFibres(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK");
    }
}

void FourChamberView::VentricularFibres(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK");
    }
}

void FourChamberView::UVCs(){
    int reply = Ask("Question", "Is this the mesh for which you would like to calculate UVCs?");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK!");
    }
    if(reply==QMessageBox::No){
        QMessageBox::information(NULL, "Answer", "Choose a different mesh.");
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
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::SplitPulmVeins(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::Corrections(){
    int reply = Ask("Question", "Does your segmentation require any of the following corrections? List options here.");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK, thanks!");
    }
}

void FourChamberView::SelectPointsCylinders() {
    if (m_Controls.button_select_pts_a->text() == "Points for Cylinders") {
        CreateInteractorWithOptions(ManualPointsType::CYLINDERS);
        m_Controls.button_select_pts_a->setText("Run Scripts for Cylinders");
    }
    else {
        // Run cylinders script
        QMessageBox::information(NULL, "Attention", "Running Scripts for Cylinders");
        m_Controls.button_select_pts_a->setEnabled(false);
    }
}

void FourChamberView::SelectPointsSlicers() {
    if (m_Controls.button_select_pts_b->text() == "Points for Slicers") {
        CreateInteractorWithOptions(ManualPointsType::SLICERS);
        m_Controls.button_select_pts_b->setText("Run Scripts for Slicers");
    }
    else {
        // Run Slicers script
        QMessageBox::information(NULL, "Attention", "Running Scripts for Slicers");
        m_Controls.button_select_pts_b->setEnabled(false);
    }
}

void FourChamberView::SelectPointsValvePlains(){
    if (m_Controls.button_select_pts_c->text() == "Points for Valve Plains") {
        CreateInteractorWithOptions(ManualPointsType::VALVE_PLAINS);
        m_Controls.button_select_pts_c->setText("Run Scripts for Valve Plains");
    }
    else {
        // Run Valve Plains script
        QMessageBox::information(NULL, "Attention", "Running Scripts for Valve Plains");
        m_Controls.button_select_pts_c->setEnabled(false);
    }
}

// helper
bool FourChamberView::RequestProjectDirectoryFromUser() {

    bool succesfulAssignment = true;

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {

        MITK_INFO << "Directory is empty. Requesting user for directory.";
        directory = QFileDialog::getExistingDirectory( NULL, "Open Project Directory",
            mitk::IOUtil::GetProgramPath().c_str(),QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

        MITK_INFO << ("Directory selected:" + directory).toStdString();

        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            MITK_WARN << "Please select a project directory with no spaces in the path!";
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            succesfulAssignment = false;
        }//_if

        if (succesfulAssignment){
            QString now = QDate::currentDate().toString(Qt::ISODate);
            QString logfilename = directory + "/LOGFILE_" + now + ".log";
            std::string logfname_str = logfilename.toStdString();

            mitk::LoggingBackend::SetLogFile(logfname_str.c_str());
            MITK_INFO << ("Changed logfile location to: " + logfilename).toStdString();
        }

    } else {
        MITK_INFO << ("Project directory already set: " + directory).toStdString();
    }//_if


    return succesfulAssignment;
}

int FourChamberView::Ask(std::string title, std::string msg){
    MITK_INFO << "[USER_INPUT]" + msg;
    return QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);
}

void FourChamberView::Warn(std::string title, std::string msg){
    MITK_WARN << msg;
    QMessageBox::warning(NULL, title.c_str(), msg.c_str());
}

QStringList FourChamberView::GetPointLabelOptions(ManualPointsType mpt) {
    QStringList res = QStringList();

    switch (mpt) {
        case ManualPointsType::CYLINDERS :
            res  << "SVC_1" << "SVC_2" << "SVC_3"  
                 << "IVC_1" << "IVC_2" << "IVC_3"  
                 << "Ao_1"  << "Ao_2"  << "Ao_3"  
                 << "PArt_1" << "PArt_2" << "PArt_3";
            break;

        case ManualPointsType::SLICERS : 
            res << "SVC_slicer_1" << "SVC_slicer_2" << "SVC_slicer_3" << "SVC_tip"
                << "IVC_slicer_1" << "IVC_slicer_2" << "IVC_slicer_3" << "IVC_tip"
                << "Ao_tip" << "PArt_tip";
            break;

        case ManualPointsType::VALVE_PLAINS :
            res << "Ao_WT_tip" << "PArt_WT_tip";
            break;
        default:
            break;
        return res;
    }

    return res;
}

void FourChamberView::CreateInteractorWithOptions(ManualPointsType mpt) {

    if (!RequestProjectDirectoryFromUser()) return;

    QStringList opts = GetPointLabelOptions(mpt);
    QString opt = GetPointTypeString(mpt);

    std::string pset_name = (opt + "_point_set").toStdString();

    // Create pointset
    mitk::PointSet::Pointer pset = mitk::PointSet::New();
    CemrgCommonUtils::AddToStorage(pset, pset_name, this->GetDataStorage(), false);

    QMessageBox::information(NULL, "Controls", "SHIFT + Left Click to add points");

    CemrgDataInteractor::Pointer m_interactor = CemrgDataInteractor::New();
    m_interactor->Initialise(opts, directory + "/" + FourChamberView::POINTS_FILE);
    m_interactor->LoadStateMachine("PointSet.xml");
    m_interactor->SetEventConfig("PointSetConfig.xml");
    m_interactor->SetDataNode(this->GetDataStorage()->GetNamedNode(pset_name.c_str()));
}

void FourChamberView::InitialiseJsonObjects() {

    QStringList pt_keys = GetPointLabelOptions(ManualPointsType::CYLINDERS);
    pt_keys.append(GetPointLabelOptions(ManualPointsType::SLICERS));
    pt_keys.append(GetPointLabelOptions(ManualPointsType::VALVE_PLAINS));

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
    fo << "#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -" << std::endl;
    fo << "#Data file for meshtools3d utility"<< std::endl;
    fo << "#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -"<< std::endl;
    fo << "[segmentation]"<< std::endl;
    fo << "seg_dir = " << meshing_parameters.seg_dir.toStdString() << std::endl;
    fo << "seg_name = " << meshing_parameters.seg_name.toStdString() << std::endl;
    fo << "mesh_from_segmentation = " << meshing_parameters.mesh_from_segmentation << std::endl;
    fo << "boundary_relabelling = " << meshing_parameters.boundary_relabelling << std::endl;
    fo << "[meshing]"<< std::endl;
    fo << "facet_angle         = "<< meshing_parameters.facet_angle << std::endl;
    fo << "facet_size          = "<< meshing_parameters.facet_size << std::endl;
    fo << "facet_distance      = "<< meshing_parameters.facet_distance << std::endl;
    fo << "cell_rad_edge_ratio = "<< meshing_parameters.cell_rad_edge_ratio << std::endl;
    fo << "cell_size           = "<< meshing_parameters.cell_size << std::endl;
    fo << "rescaleFactor       = "<< meshing_parameters.rescaleFactor << std::endl;
    fo << "[laplacesolver]" << std::endl;
    fo << "abs_toll        = " << meshing_parameters.abs_tol << std::endl;
    fo << "rel_toll        = " << meshing_parameters.rel_tol << std::endl;
    fo << "itr_max         = " << meshing_parameters.itr_max << std::endl;
    fo << "dimKrilovSp     = " << meshing_parameters.dimKrilovSp << std::endl;
    fo << "verbose         = " << meshing_parameters.verbose << std::endl;
    fo << "[others]"<< std::endl;
    fo << "eval_thickness  = " << meshing_parameters.eval_thickness << std::endl;
    fo << "[output]"<< std::endl;
    fo << "outdir          = " << meshing_parameters.out_dir.toStdString() << std::endl;
    fo << "name            = " << meshing_parameters.out_name.toStdString() << std::endl;
    fo << "out_medit       = " << meshing_parameters.out_medit << std::endl;
    fo << "out_carp        = " << meshing_parameters.out_carp << std::endl;
    fo << "out_carp_binary = " << meshing_parameters.out_carp_binary << std::endl;
    fo << "out_vtk         = " << meshing_parameters.out_vtk << std::endl;
    fo << "out_vtk_binary  = " << meshing_parameters.out_vtk_binary << std::endl;
    fo << "out_potential   = " << meshing_parameters.out_potential << std::endl;
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
            QMessageBox::information(NULL, "Attention", "No output format selected.");
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

void FourChamberView::M3dBrowseFile(const QString& dir){
    QString titlelabel, input = "";
    std::string msg;

    msg = "Select segmentation image for meshing";

    QMessageBox::information(NULL, "Attention", msg.c_str());
    input = QFileDialog::getOpenFileName(NULL, msg.c_str(), dir, QmitkIOUtil::GetFileOpenFilterString());

    QFileInfo fi(input);
    m_m3d.lineEdit_input_path->setText(input);

}
