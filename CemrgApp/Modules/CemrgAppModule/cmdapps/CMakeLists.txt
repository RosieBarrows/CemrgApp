option(BUILD_CemrgCMDApps "Build ALL command-line apps for Cemrg Module" OFF)
option(BUILD_CEMRG_FIXSHELL "Build fixing shell values command line app. " ON)
option(BUILD_CEMRG_PROJECTION_OPTS "Build Scar Map Projection Options command line app. " ON)
option(BUILD_CEMRG_CLIPPINGTOOL "Build clipping valve command line app. " ON)
option(BUILD_CEMRG_PROJECTLGE "Build lge projection command line app. " ON)
option(BUILD_CEMRG_IM2INR "Build image to inr command line app. " ON)
option(BUILD_CEMRG_VENTRICLE_SEGMENTATION_RELABEL "Build ventricle segmentation relabelling command line app" ON)
option(BUILD_CEMRG_MORPH_ANALYSIS "Build atrial morph analysis" ON)

if(BUILD_CemrgCMDApps)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgTemplateApp
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgAppCmdTemplate.cpp
  )
endif()

if(BUILD_CEMRG_FIXSHELL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgFixShell
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgFixShellApp.cpp
  )
endif()

if(BUILD_CEMRG_PROJECTION_OPTS)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgScarProjectionOptions
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgProjectionOptions.cpp
  )
endif()

if(BUILD_CEMRG_CLIPPINGTOOL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgClippingTool
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgClippingTool.cpp
  )
endif()

if(BUILD_CEMRG_PROJECTLGE)
  mitkFunctionCreateCommandLineApp(
    NAME ProjectLgeToVtkMesh
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgProjectLgeToVtkMesh.cpp
  )
endif()

if(BUILD_CEMRG_IM2INR)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgIM2INR
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgIM2INR.cpp
  )
endif()

if(BUILD_CEMRG_VENTRICLE_SEGMENTATION_RELABEL)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgVentricleSegmRelabel
    DEPENDS MitkCemrgAppModule
    CPP_FILES CemrgVentricleSegmRelabel.cpp
  )
endif()

if(BUILD_CEMRG_MORPH_ANALYSIS)
  mitkFunctionCreateCommandLineApp(
    NAME CemrgMorphAnalysis
    DEPENDS MitkSegmentation MitkCemrgAppModule
    CPP_FILES CemrgMorphAnalysis.cpp
  )
endif()
