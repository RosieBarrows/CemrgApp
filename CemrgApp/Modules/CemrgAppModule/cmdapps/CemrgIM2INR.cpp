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
CEMRG CMD APP TEMPLATE
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkSurface.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkImplicitBoolean.h>
#include <vtkImplicitVolume.h>
#include <vtkImplicitDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkDataSetSurfaceFilter.h>

// ITK
#include <itkPoint.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageFileWriter.h>

// Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>

// C++ Standard
#include <algorithm>
#include <string>
#include <numeric>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgCommonUtils.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("Image to INR Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription(
        "Convert an image file (e.g NII) to an INR.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::InputFile,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
        "input", "i", mitkCommandLineParser::InputFile,
        "Input image", "Full path of image file.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output file", "Where to save the output (e.g output.inr).");
    parser.addArgument(
        "uint8", "u", mitkCommandLineParser::Bool,
        "Convert to UINT8", "Convert image type to UINT8 (default=false).");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    // auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    auto verbose = false;
    auto convert2uint = false;
    std::string outFilename = "convert.inr";

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose"))
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);

    if (parsedArgs.end() != parsedArgs.find("uint8"))
        convert2uint = us::any_cast<bool>(parsedArgs["uint8"]);

    if (parsedArgs.end() != parsedArgs.find("output"))
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);


    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(inFilename);
        QString outname = QString::fromStdString(outFilename);

        if (!outname.contains(".inr", Qt::CaseSensitive))
            outname = outname + ".inr";

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and inputPath variables
        QFileInfo fi(inname);

        QString direct = fi.absolutePath();
        QString inputPath = fi.absoluteFilePath();
        // QString outputPath = direct + "/" + outname;

        QString outpath = CemrgCommonUtils::ConvertToInr(direct, fi.fileName(), convert2uint, outname);

        if(outpath.isEmpty()){
            MITK_ERROR << "Problems creating the file";
            return EXIT_FAILURE;
        }

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
