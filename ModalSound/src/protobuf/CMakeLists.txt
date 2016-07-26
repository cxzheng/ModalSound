# find google protobuf library
find_package(Protobuf REQUIRED)

include_directories(
    ${EIGEN_INCLUDE_DIR}
    ${PROTOBUF_INCLUDE_DIRS})
PROTOBUF_GENERATE_CPP(PROTO_SRCS PROTO_HDRS sploosh.proto)

# -----------------------------------------------------------------------------
# generate python protobuf file
set(PY_SRC "${BIN_DIR}/python/sploosh_pb2.py")
get_filename_component(ABS_PROTO_FILE sploosh.proto ABSOLUTE)
get_filename_component(ABS_PROTO_PATH ${ABS_PROTO_FILE} PATH)
add_custom_command(
    OUTPUT ${PY_SRC}
    COMMAND  ${PROTOBUF_PROTOC_EXECUTABLE}
    ARGS --python_out ${BIN_DIR}/python -I ${ABS_PROTO_PATH} ${ABS_PROTO_FILE}
    DEPENDS ${ABS_PROTO_FILE}
    COMMENT "Running python protocol buffer compiler on ${ABS_PROTO_FILE}"
    VERBATIM)

add_exe(gen_moments
    gen_moments.cpp
    ${PROTO_HDRS}   ${PROTO_SRCS}
    ../generic/precision_type.hpp
    ../generic/precision_type.cpp
    ${PY_SRC}
    LINK_LIBRARIES ${Boost_SYSTEM_LIBRARY}
                   ${Boost_FILESYSTEM_LIBRARY}
                   ${Boost_PROGRAM_OPTIONS_LIBRARY}
                   ${GSL_LIBRARIES}
                   ${PROTOBUF_LIBRARY}
    OUT_DIR ${BIN_DIR})
