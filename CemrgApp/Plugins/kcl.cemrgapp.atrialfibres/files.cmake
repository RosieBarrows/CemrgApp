set(SRC_CPP_FILES

)

set(INTERNAL_CPP_FILES
  kcl_cemrgapp_atrialfibres_Activator.cpp
  AtrialFibresView.cpp
  AtrialFibresClipperView.cpp
  AtrialFibresLandmarksView.cpp
  AtrialFibresVisualiseView.cpp
)

set(UI_FILES
  src/internal/AtrialFibresViewControls.ui
  src/internal/AtrialFibresViewUIAnalysisSelector.ui
  src/internal/AtrialFibresViewUIEditLabels.ui
  src/internal/AtrialFibresViewUIUacSelector.ui
  src/internal/AtrialFibresViewUIMeshing.ui
  src/internal/AtrialFibresViewUIRemesh.ui
  src/internal/AtrialFibresViewUIConvert.ui
  ../kcl.cemrgapp.scar/src/internal/AtrialScarViewUIcemrgnet.ui
  src/internal/AtrialFibresClipperViewControls.ui
  src/internal/AtrialFibresClipperViewLabels.ui
  src/internal/AtrialFibresClipperViewUIRadius.ui
  src/internal/AtrialFibresClipperViewUICorridor.ui
  src/internal/AtrialFibresLandmarksViewControls.ui
  src/internal/AtrialFibresLandmarksViewRough.ui
  src/internal/AtrialFibresLandmarksViewRefined.ui
  src/internal/AtrialFibresVisualiseViewControls.ui
)

set(MOC_H_FILES
  src/internal/kcl_cemrgapp_atrialfibres_Activator.h
  src/internal/AtrialFibresView.h
  src/internal/AtrialFibresClipperView.h
  src/internal/AtrialFibresLandmarksView.h
  src/internal/AtrialFibresVisualiseView.h
)

# list of resource files which can be used by the plug-in
# system without loading the plug-ins shared library,
# for example the icon used in the menu and tabs for the
# plug-in views in the workbench
set(CACHED_RESOURCE_FILES
  resources/icon.xpm
  plugin.xml
)

# list of Qt .qrc files which contain additional resources
# specific to this plugin
set(QRC_FILES

)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
