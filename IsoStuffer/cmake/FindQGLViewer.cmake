IF (QGLVIEWER_ROOT)                                                                                               
    FIND_PATH(QGLVIEWER_INCLUDE_DIR QGLViewer/qglviewer.h                                                         
        HINTS                                                                                                     
            PATHS "/usr/local/include"
            PATHS ${QGLVIEWER_ROOT}/include)                                                                      
    FIND_LIBRARY(QGLVIEWER_LIBRARY QGLViewer                                                                      
        HINTS                                                                                                     
            PATHS ${QGLVIEWER_ROOT}                                                                               
            PATHS ${QGLVIEWER_ROOT}/lib                                                                           
            PATHS ${QGLVIEWER_ROOT}/Release)                                                                      
ELSE (QGLVIEWER_ROOT)                                                                                             
  FIND_PATH(QGLVIEWER_INCLUDE_DIR qglviewer.h                                                                   
        HINTS                                                                                                     
            PATHS "/usr/local/include/QGLViewer"
            PATHS "/usr/include/QGLViewer"
            PATHS "/usr/local/lib/QGLViewer.framework/Versions/Current/Headers/"
            PATHS "/usr/include/QGLViewer/"
            PATHS ${SYSTEM_INC_PATH}                             
            PATHS $ENV{INCLUDE}
            )
    FIND_LIBRARY(QGLVIEWER_LIBRARY QGLViewer-qt5                                                                     
      HINTS
            PATHS "/usr/local/lib/"
            PATHS "/usr/lib/"
            PATHS "/usr/local/lib/QGLViewer.framework/Versions/Current/"
            PATHS ${SYSTEM_LIB_PATH}                                                                              
            PATHS $ENV{LD_LIBRARY_PATH}
            )
    if (QGLVIEWER_LIBRARY)
      message("QGLViewer library: ${QGLVIEWER_LIBRARY}")
    else (QGLVIEWER_LIBRARY)
      FIND_LIBRARY(QGLVIEWER_LIBRARY QGLViewer 
        HINTS
              PATHS "/usr/local/lib/"
              PATHS "/usr/lib/"
              PATHS "/usr/local/lib/QGLViewer.framework/Versions/Current/"
              PATHS ${SYSTEM_LIB_PATH}                                                                       
              PATHS $ENV{LD_LIBRARY_PATH}
              )
          endif(QGLVIEWER_LIBRARY)


ENDIF (QGLVIEWER_ROOT)                                                                                            
                                                                                                                  
MARK_AS_ADVANCED(QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY)                                                         
message("QGLViewer include: ${QGLVIEWER_INCLUDE_DIR}")
                                                                                                                  
INCLUDE(FindPackageHandleStandardArgs)                                                                            
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QGLViewer DEFAULT_MSG                                                           
    QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY)  
