# - Find CPPUnit
FIND_PATH(CppUnit_DIR cppunit/Test.h
    PATHS ${SYSTEM_INC_PATH}
    PATHS $ENV{INCLUDE})

FIND_LIBRARY(CppUnit_LIB cppunit
    PATHS ${SYSTEM_LIB_PATH}
    PATHS $ENV{LD_LIBRARY_PATH})

IF (CppUnit_DIR AND CppUnit_LIB)
    IF (NOT CppUnit_FIND_QUIETLY)
        MESSAGE(STATUS "Found CppUnit: ${CppUnit_LIB}")
    ENDIF (NOT CppUnit_FIND_QUIETLY)
ELSE (CppUnit_DIR AND CppUnit_LIB)
    IF (CppUnit_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find CppUnit")
    ENDIF (CppUnit_FIND_REQUIRED)
ENDIF (CppUnit_DIR AND CppUnit_LIB)

