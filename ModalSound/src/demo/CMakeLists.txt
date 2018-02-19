# find google protobuf library
find_package(Protobuf REQUIRED)

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

PROTOBUF_GENERATE_CPP(PROTO_SRCS PROTO_HDRS ../protobuf/sploosh.proto)

include_directories(
  ${CMAKE_INCLUDE_CURRENT_DIR}
  )

message("QGLGLGLGLG: ${QGLVIEWER_LIBRARY}")

add_exe(click_synth
    ClickSoundSynth.cpp
    AudioProducer.h     AudioProducer.cpp
    ../modal/ModalModel.h
    ../modal/ModalModel.cpp
    ../transfer/FMMTransferEval.h
    ../transfer/FMMTransferEval.cpp
    ../io/MatrixIO.hpp
    ../io/MatrixIO.cpp
    ../generic/precision_type.hpp
    ../generic/precision_type.cpp
    ${PROTO_HDRS}   ${PROTO_SRCS}
    LINK_LIBRARIES  ${GLUT_LIBRARIES}   
                    ${OPENGL_LIBRARIES}
                    ${GLUT_glut_LIBRARY}
                    ${GLUT_Xi_LIBRARY}
                    ${GLUT_Xmu_LIBRARY}
                    ${OPENGL_glu_LIBRARY}
                    ${OPENGL_gl_LIBRARY}
                    Qt5::Core   Qt5::Widgets
                    Qt5::Xml    Qt5::OpenGL
                    Qt5::Multimedia
                    ${MKL_LIBS}
                    ${QGLVIEWER_LIBRARY}
                    ${PROTOBUF_LIBRARY}
                    ${Boost_FILESYSTEM_LIBRARY} 
                    ${Boost_SYSTEM_LIBRARY}
                    ${GSL_LIBRARIES}
    OUT_DIR ${BIN_DIR}
)
