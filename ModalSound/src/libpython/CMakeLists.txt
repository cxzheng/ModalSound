include_directories(
    ${Eigen_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR})

# source files
set(SRCS
    export.cpp
    ../modal/ModalShape.h
    ../modal/ModalShape.cpp
    ../generic/precision_type.hpp
    ../generic/precision_type.cpp
)

add_library(sploosh MODULE ${SRCS})

target_link_libraries(sploosh
    ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${Boost_PYTHON_LIBRARY}     ${PYTHON_LIBRARY}
    ${GSL_LIBRARIES}
)

set_target_properties(sploosh PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR}/python)

set_target_properties(sploosh PROPERTIES OUTPUT_NAME sploosh PREFIX "")
if (APPLE)
    # Python seems to only look for plugins that end in .so, not .dylib                                     
    set_target_properties(sploosh PROPERTIES SUFFIX ".so")
elseif (WIN32)
    # On Windows, Python looks for .pyd instead of .dll                                                     
    set_target_properties(sploosh PROPERTIES SUFFIX ".pyd")
endif()
