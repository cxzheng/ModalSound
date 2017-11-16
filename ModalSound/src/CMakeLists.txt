# utils/
set(UTILS_SRC
    utils/macros.h      utils/math.hpp)
# linearalgebra
set(SC_SRC
    sc/Quaternion.hpp   sc/Quaternion.cpp)

include_directories(
    ${Eigen_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${MKL_INCLUDE_DIR})

###############################################################################
if ( BUILD_TOOLS )
    add_subdirectory(tools)
endif ( BUILD_TOOLS )

if ( BUILD_PYTHON )
    add_subdirectory(libpython)
endif ( BUILD_PYTHON )

if ( BUILD_DEMO )
    add_subdirectory(demo)
endif ( BUILD_DEMO )

add_subdirectory(protobuf)
