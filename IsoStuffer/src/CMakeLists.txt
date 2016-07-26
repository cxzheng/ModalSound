## CMake specification for src directory
set(SOURCES config.h
    tables.h            IsoStuffer.hpp
    OctTree.h           OctTree.cpp
    DistValProc.h       DistValProc.cpp)

# utils/
set(UTILS_SRC
    utils/macros.h      utils/math.hpp)

# linearalgebra
set(LINEAR_SRC
    linearalgebra/eig3.h    linearalgebra/eig3.cpp)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Qt5Core_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5Xml_INCLUDE_DIRS}
    ${Qt5OpenGL_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

if ( USE_GUI )
    include_directories(
        ${GLUT_INCLUDE_DIR}
        ${OPENGL_INCLUDE_DIR})

    if ( NOT QGLVIEWER_FRAMEWORK )
        include_directories(${QGLVIEWER_INCLUDE_DIR})
    endif ()

    # Find includes in corresponding build directories
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTORCC ON)

    # ====================================================================

    # setup the qrc files
    set(UI_QRC_SRCS ui/images.qrc)
    ##QT5_ADD_RESOURCES(UI_QRC_SRCS ${UI_QRC_FILES})

    # headers to be moc-ified
    set(UI_MOC_HEADERS 
        ui/IsoStufferFrame.h    ui/IsoStufferCanvas.h
        ui/ParamsDialog.h)
    # the ui files
    set(UI_UI_FILES ui/IsoStufferFrame.ui ui/Params.ui)

    QT5_WRAP_CPP(UI_MOC_SRCS ${UI_MOC_HEADERS})
    ##    OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    QT5_WRAP_UI(UI_UI_SRCS ${UI_UI_FILES})

    set(GUIAPP_SRCS GuiMain.cpp
        ui/IsoStufferFrame.h    ui/IsoStufferFrame.cpp
        io/TglMeshReader.hpp    io/TglMeshReader.cpp
        ui/IsoStufferCanvas.h   ui/IsoStufferCanvas.cpp
        generic/precision_type.hpp
        generic/precision_type.cpp
        ${UI_UI_SRCS} ${UI_MOC_SRCS} 
        ${UI_QRC_SRCS} )
    
    # =========================================================================
    # set for the final executable for isostuffer
    add_executable(isostuffer-gui
        ${GUIAPP_SRCS}          ${SOURCES}              ${UTILS_SRC}
        ${LINEAR_SRC})

    target_link_libraries(isostuffer-gui
        ${GLUT_LIBRARIES}       
        Qt5::Core   Qt5::Widgets
        Qt5::Xml    Qt5::OpenGL
        ${QGLVIEWER_LIBRARY}
        ${OPENGL_LIBRARIES}     
        ${OPENMP_LIBS}          
        ${ICC_LIBS}
    )

    set_target_properties(isostuffer-gui PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( USE_GUI ) 


# =============================================================================
add_executable(isostuffer
    io/TglMeshReader.hpp    io/TglMeshReader.cpp
    generic/precision_type.hpp
    generic/precision_type.cpp
    main.cpp
    ${SOURCES}              ${UTILS_SRC}            ${LINEAR_SRC})

target_link_libraries(isostuffer ${ICC_LIBS}        
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

set_target_properties(isostuffer PROPERTIES LINK_FLAGS ${LINK_FLAGS})

# # =============================================================================
# add_executable(ellip_gen
#     ellip_gen.cpp               ${SOURCES}              ${UTILS_SRC}
#     generic/precision_type.cpp
#     ${LINEAR_SRC})
# target_link_libraries(ellip_gen ${ICC_LIBS})
# set_target_properties(ellip_gen PROPERTIES LINK_FLAGS ${LINK_FLAGS})
# 
# add_executable(ellip_gen_abc
#     ellip_gen_abc.cpp           ${SOURCES}              ${UTILS_SRC}
#     generic/precision_type.cpp
#     ${LINEAR_SRC})
# target_link_libraries(ellip_gen_abc ${ICC_LIBS} ${OPENMP_LIBS})
# set_target_properties(ellip_gen_abc PROPERTIES LINK_FLAGS ${LINK_FLAGS})

