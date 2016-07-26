# - Find library libConfig++ 
IF ( LibConfig_DIR )
    FIND_PATH(LibConfig_INCLUDE_DIR libconfig.h++
        PATHS ${LibConfig_DIR}/include)
    FIND_LIBRARY(LibConfig_LIBRARY config++
        PATHS ${LibConfig_DIR}/lib)
ELSE ( LibConfig_DIR )
    FIND_PATH(LibConfig_INCLUDE_DIR libconfig.h++
        PATHS ${SYSTEM_INC_PATH}
        PATHS $ENV{INCLUDE})
    FIND_LIBRARY(LibConfig_LIBRARY config++
        PATHS ${SYSTEM_LIB_PATH}
        PATHS $ENV{INCLUDE})
ENDIF ( LibConfig_DIR )

IF (LibConfig_INCLUDE_DIR AND LibConfig_LIBRARY)
    IF (NOT LibConfig_FIND_QUIETLY)
        MESSAGE(STATUS "Found libConfig: ${LibConfig_LIBRARY}")
    ENDIF (NOT LibConfig_FIND_QUIETLY)
ELSE (LibConfig_INCLUDE_DIR AND LibConfig_LIBRARY)
    IF (LibConfig_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find libConfig")
    ENDIF(LibConfig_FIND_REQUIRED)
ENDIF (LibConfig_INCLUDE_DIR AND LibConfig_LIBRARY)

