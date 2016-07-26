find_package(Boost 1.50 OPTIONAL_COMPONENTS "python") 

# Python libraries
find_package(PythonLibs "2.7")
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON 
    "Build python bindings." ON
    "PYTHONLIBS_FOUND;Boost_PYTHON_LIBRARY" OFF)

