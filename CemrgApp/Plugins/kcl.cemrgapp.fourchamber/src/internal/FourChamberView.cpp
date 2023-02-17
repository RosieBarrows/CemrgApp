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
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkImagePixelReadAccessor.h>

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
#include <QDate>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>


const std::string FourChamberView::VIEW_ID = "org.mitk.views.fourchamberheart";

void FourChamberView::SetFocus(){
    // m_Controls.buttonPerformImageProcessing->setFocus();
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
    connect(m_Controls.buttonPerformImageProcessing, SIGNAL(clicked()), this, SLOT(SetFocus()));

    connect(m_Controls.button_select_pts_a, SIGNAL(clicked()), this, SLOT(SelectPointsA()));
    connect(m_Controls.button_select_pts_b, SIGNAL(clicked()), this, SLOT(SelectPointsB()));

    // Set default variables and initialise objects
    m_Controls.button_loaddicom->setVisible(false);
    m_Controls.button_extractsurfs->setVisible(false);
    m_Controls.button_uvclandmarks->setVisible(false);
    m_Controls.button_calcuvcs->setVisible(false);
}

void FourChamberView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

// SLOTS
void FourChamberView::SetWorkingFolder(){
    if (!RequestProjectDirectoryFromUser()){
        MITK_WARN << "Folder not set. Check LOGFILE."; 
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
            MITK_INFO << ("Conversion succesful. Intermediate NII folder: "  tmpNiftiFolder).toStdString();
            QMessageBox::information(NULL, "Information", "Conversion successful, press the Process Images button to continue.");
            QDir niftiFolder(tmpNiftiFolder);
            QStringList niftiFiles = niftiFolder.entryList();

            if (niftiFiles.size()>0) {

                QString thisFile, path;
                for(int ix=0; ix<niftiFiles.size(); ix+) {

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

void FourChamberView::Meshing(){

    if (!RequestProjectDirectoryFromUser()) return; 

    QString path = FourChamberView::directory + "/segmentations";
    QString segname = "seg_final_smooth.nrrd";
    path += "/" + segname;

    if(QFile::exists(path)){
        std::string msg = "Load previously smoothed segmentation (" + segname.toStdString() + ")?";
        int reply = QMessageBox::question(NULL, "Question", msg.c_str());
        if(reply==QMessageBox::Yes){
            MITK_INFO << "Loading segmentation from file";
            // do stuff! 
        } else {
            // retrieve the path of the segmentation after user has selected it from the UI
            path = QFileDialog::getOpenFileName(NULL, "Open Segmentation file",
                directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        }
    }

    QString parpath = QFileDialog::getOpenFileName(NULL, "Open Parameters file (.par)",
    directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

    // convert segmentation to inr
    try {
        // load segmentation 
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
        int dimensions = image->GetDimension(0) * image->GetDimension(1) * image->GetDimension(2);

        //Convert image to right type
        // itk::Image<uint8_t, 3>::Pointer itkImage = itk::Image<uint8_t, 3>::New();
        // mitk::CastToItkImage(image, itkImage);
        // mitk::CastToMitkImage(itkImage, image);

        //Access image volume
        mitk::ImagePixelReadAccessor<uint8_t, 3> readAccess(image);
        uint8_t* pv = (uint8_t*)readAccess.GetData();

        //Prepare header of inr file (BUGS IN RELEASE MODE DUE TO NULL TERMINATOR \0)
        char header[256] = {};
        int bitlength = 8;
        const char* btype = "unsigned fixed";
        mitk::Vector3D spacing = image->GetGeometry()->GetSpacing();
        int n = sprintf(header, "#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\nPIXSIZE=%d bits\nCPU=decm\nVX=%6.4f\nVY=%6.4f\nVZ=%6.4f\n", image->GetDimension(0), image->GetDimension(1), image->GetDimension(2), btype, bitlength, spacing.GetElement(0), spacing.GetElement(1), spacing.GetElement(2));
        for (int i = n; i < 252; i++)
            header[i] = '\n';

        header[252] = '#';
        header[253] = '#';
        header[254] = '}';
        header[255] = '\n';

        //Write to binary file
        QString inr_path = directory + "/segmentations/converted.inr";
        std::string path = inr_path.toStdString();
        ofstream myFile(path, ios::out | ios::binary);
        myFile.write((char*)header, 256 * sizeof(char));
        myFile.write((char*)pv, dimensions * sizeof(uint8_t));
        myFile.close();

        // create cmd object (CemrgCommandLine) outputs to directory/meshing
        std::unique_ptr<CemrgCommandLine> cmd_object(new CemrgCommandLine());
        QString meshing_output = cmd_object->ExecuteCreateCGALMesh(directory, "meshname", parpath, inr_path, "meshing");

    } catch (mitk::Exception& e) {
        //Deal with the situation not to have access
        qDebug() << e.GetDescription();
        return;
    }//_try
    
    MITK_INFO << path.toStdString();


    // int reply = Ask("Question", "Are these the correct labels for your segmentation?");
    // if(reply==QMessageBox::Yes){
    //     QMessageBox::information(NULL, "Answer", "OK!");
    // }
    // if(reply==QMessageBox::No){
    //     QMessageBox::information(NULL, "Answer", "Manually change labels from default OR select button to provide a .json file");
    // }   
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

void FourChamberView::PrepareSegmentation(){
    int reply = Ask("Question", "Is this the segmentation you would like to use to create your heart model?");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::CalculateUVCs(){
    int reply = Ask("Question", "Placeholder");
    if(reply==QMessageBox::Yes){
        QMessageBox::information(NULL, "Answer", "OK!");
    }
}

void FourChamberView::SelectPointsA() {
    CreateInteractorWithOptions("A");
}

void FourChamberView::SelectPointsB() {
    CreateInteractorWithOptions("B");
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
    return QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);
}

QStringList FourChamberView::GetPointLabelOptions(QString opt) {
    QStringList res = QStringList();
    if (opt == "A") {
        res  << "SVC_1"  
             << "SVC_2"  
             << "SVC_3"  
             << "IVC_1"  
             << "IVC_2"  
             << "IVC_3"  
             << "Ao_1"  
             << "Ao_2"  
             << "Ao_3"  
             << "PArt_1"  
             << "PArt_2"  
             << "PArt_3";
    }
    else if (opt == "B") {
        res << "SVC_slicer_1"
            << "SVC_slicer_2"
            << "SVC_slicer_3"
            << "SVC_tip"
            << "IVC_slicer_1"
            << "IVC_slicer_2"
            << "IVC_slicer_3"
            << "IVC_tip"
            << "Ao_tip"
            << "PArt_tip";
    }else if (opt == "C") {
        res << "Ao_WT_tip"
            << "PArt_WT_tip";
    }
    return res;
}

void FourChamberView::CreateInteractorWithOptions(QString opt) {
    QStringList opts = GetPointLabelOptions(opt);
      // Create pointset
      mitk::PointSet::Pointer pset = mitk::PointSet::New();
    pset->SetName(opt  "_point_set");
      // create node
      mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(pset);
    node->SetName(opt  "_point_set");
    node->SetVisibility(true);
     CemrgDataInteractor m_interactor = CemrgDataInteractor::New(opts);
    m_interactor->LoadStateMachine("PointSet.xml");
    m_interactor->SetEventConfig("PointSetConfig.xml");
    m_interactor->SetDataNode(node);
    
}