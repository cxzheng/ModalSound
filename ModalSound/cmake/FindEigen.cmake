find_path(Eigen_INCLUDE_DIR Eigen/Core                          
  HINTS
     PATHS "/usr/local/include/eigen3"
     PATHS "/usr/include/eigen3"
     PATHS "/usr/local/Cellar/eigen/3.2.4/include/eigen3"
     PATHS "/usr/local/Cellar/eigen/3.2.6/include/eigen3"
     PATHS ${EIGEN_ROOT}                                           
     PATHS "~/Codes/eigen/")

mark_as_advanced(Eigen_INCLUDE_DIR)                                  
include(FindPackageHandleStandardArgs)                               
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen DEFAULT_MSG Eigen_INCLUDE_DIR) 
