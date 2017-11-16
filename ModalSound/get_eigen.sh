#!/bin/bash

actual_eigen_tar_md5="a7aab9f758249b86c93221ad417fbe18"
eigen_url="http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2"
eigen_file_name="3.3.4.tar.bz2"

extracted_eigen_name="eigen-eigen-5a0156e40feb"
# md5 on installed Eigen source files
actual_installed_eigen_md5="1829213ea6385ea2a942258d91bf7877"

# Verify that curl is installed
command -v curl >/dev/null 2>&1 || { echo >&2 "Error, please install curl and rerun get_eigen.sh."; exit 1; }

# Verify that md5sum is installed
command -v md5sum >/dev/null 2>&1 || command -v md5 >/dev/null 2>&1 || { echo >&2 "Error, please install md5 or md5sum and rerun get_dependencies.sh."; exit 1; }
if command -v md5 > /dev/null 2>&1 ; then
  md5sum()
  {
    md5 "$@" | sed -e 's#^MD5 [(]\(.*\)[)] = \(.*\)$#\2 \1#'
  }
  export -f md5sum
fi

# If the output directory exists
if [ -d "include/eigen" ]; then
  # If the Eigen install is up to date, no action is needed
  computed_installed_eigen_md5=`find include/eigen -type f -name '*.h' -exec bash -c 'md5sum "$0" "$@"' {} + | awk '{print $2$1}' | sort -fd | md5sum | cut -c -32`
  if [ "$computed_installed_eigen_md5" == "$actual_installed_eigen_md5" ]
  then
    echo "Eigen matrix library is already up to date, no further action is needed."
    exit 0
  fi
  # Otherwise, the checksum is incorrect, warn the user and exit
  echo "Error, directory include/eigen has an incorrect checksum, please run remove_eigen.sh and rerun get_eigen.sh."
  exit 1
fi

echo "Installing Eigen"

# Create a temporary working directory
temp_dir_name=`uuidgen`
if [ -d "$temp_dir_name" ]; then
  echo "Error, temporary working directory $temp_dir_name exists, this is a bug. Please contact the maintainer."
  exit 1
fi
echo "--->  Creating temporary directory $temp_dir_name"
mkdir $temp_dir_name

function cleanup {
  echo "--->  Removing temporary directory $temp_dir_name"
  rm -fr "$temp_dir_name"
}
trap cleanup EXIT

# Download Eigen
echo "--->  Downloading Eigen source"
curl -s -L -o "$temp_dir_name/$eigen_file_name" "$eigen_url"
if [ $? -ne 0 ]
then
  echo "Error, failed to download Eigen from $eigen_url."
  exit 1
fi

# Run a checksum on the download
echo "--->  Verifying Eigen checksum"
computed_eigen_tar_md5=`md5sum $temp_dir_name/$eigen_file_name | cut -c -32`
if [ "$actual_eigen_tar_md5" != "$computed_eigen_tar_md5" ]
then
  echo "Error, md5 checksum for $eigen_file_name does not match $actual_eigen_tar_md5."
  exit 1
fi

# Extract the tar archive
echo "--->  Extracting Eigen"
tar -xf "$temp_dir_name"/"$eigen_file_name" -C "$temp_dir_name"
# Move the source to its final location
echo "--->  Moving Eigen to destination"
mkdir -p include/eigen
mv $temp_dir_name/$extracted_eigen_name/Eigen include/eigen/
if [ $? -ne 0 ]
then
  echo "Error, failed to move Eigen source from $temp_dir_name/$extracted_eigen_name/Eigen to include/eigen/."
  exit 1
fi
mv $temp_dir_name/$extracted_eigen_name/signature_of_eigen3_matrix_library include/eigen/
if [ $? -ne 0 ]
then
  echo "Error, failed to move Eigen signature from $temp_dir_name/$extracted_eigen_name/signature_of_eigen3_matrix_library to include/eigen/."
  exit 1
fi

trap - EXIT
cleanup
echo "Successfully installed Eigen"
