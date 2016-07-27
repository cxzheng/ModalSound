include_directories(
    ${GLUT_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR})

if ( NOT QGLVIEWER_FRAMEWORK )
    include_directories(${QGLVIEWER_INCLUDE_DIR})
endif ()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTORCC ON)

# =================== TETVIEWER ==================
# the qrc files
set(TETVIEWER_QRC_SRCS images.qrc)

QT5_WRAP_UI(TETVIEWER_UI_SRC 
        tetviewer.ui 
        ModesVisParams.ui)

QT5_WRAP_CPP(TETVIEWER_MOC_SRC 
        TetViewerFrame.h TetViewerCanvas.h
        ModesVisParamsDialog.h)

add_executable(tetviewer
    TetViewerFrame.h      TetViewerFrame.cpp    
    TetViewerCanvas.h     TetViewerCanvas.cpp
    ModesVisParamsDialog.h
    ../generic/precision_type.hpp  ../generic/precision_type.cpp
    ${TETVIEWER_QRC_SRCS}
    ${TETVIEWER_UI_SRC}         ${TETVIEWER_MOC_SRC})

target_link_libraries(tetviewer
    ${QGLVIEWER_LIBRARY}
    ${GLUT_LIBRARIES}           
    ${OPENGL_LIBRARIES}
    Qt5::Core   Qt5::Widgets
    Qt5::Xml    Qt5::OpenGL
)

set_target_properties(tetviewer 
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# -----------------------------------------------------------------------------
add_executable(extmat
    extract_matrices.cpp
    ../sc/PardisoMatrix.hpp        ../sc/PardisoMatrix.cpp
    ../deformable/stvk.h           ../deformable/stvk.cpp
    ../generic/precision_type.hpp  ../generic/precision_type.cpp)
target_link_libraries(extmat
    ${MKL_LIBS} ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(extmat 
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
