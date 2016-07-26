# - Find library libQGLViewer
IF (QGLVIEWER_ROOT)
    # -------------------------------------------------------
    FIND_LIBRARY(QGLVIEWER_LIBRARY QGLViewer
        HINTS
            PATHS ${QGLVIEWER_ROOT}
            PATHS ${QGLVIEWER_ROOT}/lib
            PATHS ${QGLVIEWER_ROOT}/Release)
    # -------------------------------------------------------
    FIND_PATH(QGLVIEWER_INCLUDE_DIR QGLViewer/qglviewer.h
        HINTS
            PATHS ${QGLVIEWER_ROOT}/include)

ELSE (QGLVIEWER_ROOT)

    # -------------------------------------------------------
    FIND_LIBRARY(QGLVIEWER_LIBRARY QGLViewer
        HINTS
            PATHS ${SYSTEM_LIB_PATH}
            PATHS $ENV{LD_LIBRARY_PATH})
    # -------------------------------------------------------
    FIND_PATH(QGLVIEWER_INCLUDE_DIR QGLViewer/qglviewer.h
        HINTS
            PATHS ${SYSTEM_INC_PATH}
            PATHS $ENV{INCLUDE})

ENDIF (QGLVIEWER_ROOT)

MARK_AS_ADVANCED(QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY)

INCLUDE(FindPackageHandleStandardArgs)

# ------------------------------------------------------------
# ------------------------------------------------------------

if ( APPLE AND QGLVIEWER_LIBRARY MATCHES ".*framework$" ) 

    FIND_PACKAGE_HANDLE_STANDARD_ARGS(QGLViewer DEFAULT_MSG QGLVIEWER_LIBRARY)
    set(QGLVIEWER_FRAMEWORK ON)

else ()

    FIND_PACKAGE_HANDLE_STANDARD_ARGS(QGLViewer DEFAULT_MSG 
        QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY)
    set(QGLVIEWER_FRAMEWORK OFF)

endif ()
