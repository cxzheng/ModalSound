# turn on the features of Qt for libQGLViewer
set(QT_USE_QTXML 1)
set(QT_USE_QTOPENGL 1)

include(${QT_USE_FILE})
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
    ${GLUT_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    ${QGLVIEWER_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR})

# -----------------------------------------------------------------------------
add_executable(collisionviewer
    CollisionViewer.cpp       ../sc/eig3.cpp
    ../collision/RigidCollisionDetector.cpp
    ../generic/precision_type.cpp
    ../geometry/tritri.cpp)
target_link_libraries(collisionviewer
    ${QT_LIBRARIES}             ${QGLVIEWER_LIBRARY}
    ${GLUT_LIBRARIES}           ${OPENGL_LIBRARIES})
set_target_properties(collisionviewer 
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
