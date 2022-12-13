#!/bin/bash

#-----------------------------------------------------------------------------
# This script configures, compiles and install packages bundled to UltraCold
#-----------------------------------------------------------------------------

#-------------------------------------------
# Configure, build and install arpack-ng
#-----------------------------------------------------------------------------

SUCCESS_FILE="arpack-ng-success.txt"

if [ -f "$SUCCESS_FILE" ]; then

  # Check if success file exists

  echo "-- Found the file arpack-ng/$SUCCESS_FILE, so assume that arpack-ng is already correctly configured, built and
  installed from a previous run of cmake."

else

  # Enter arpack-ng folder and create the build directory

  cd arpack-ng

  if ! [[ -d "build" ]]; then
  
    mkdir build

  fi

  # Run cmake from the CMakeLists.txt distributed with arpack-ng

  echo "-- Configuring arpack-ng"
  
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=$1 \
        -DCMAKE_Fortran_COMPILER=ifort \
        -DCMAKE_C_COMPILER=icc \
        -DCMAKE_CXX_COMPILER=icpc \
        -DICB=ON \
        -DCMAKE_BUILD_TYPE=Release \
        .. 

  # Check if make is able to build arpack-ng

  echo "-- Compiling arpack-ng..."

  fail_to_run_make=$( make -j 4 2>&1 > /dev/null )
  if ! [[ -z $fail_to_run_make ]];then

    echo -e " -- ${RED} ERROR: Failed to run make for arpack-ng ${NC}"
    echo -e " -- ${RED} The error is: $fail_to_run_make ${NC}"

    if [[ -d "$SUCCESS_FILE" ]]; then

	    rm $SUCCESS_FILE

    fi 

    exit 1 

  fi

  # Check if make is able to test arpack-ng

  echo "-- Testing arpack-ng..."

  fail_to_run_make_test=$( make test 2>&1 > /dev/null )
  if ! [[ -z $fail_to_run_make_test ]];then

    echo -e " -- ${RED} ERROR: Failed to run make test for arpack-ng ${NC}"
    echo -e " -- ${RED} The error is: $fail_to_run_make_test ${NC}" 

    if [[ -d "$SUCCESS_FILE" ]]; then

	    rm $SUCCESS_FILE

    fi 

    exit 1 

  fi

  # Check if make is able to install arpack-ng

  echo "-- Installing arpack-ng..."

  fail_to_run_make_install=$( make install 2>&1 > /dev/null )
  if ! [[ -z $fail_to_run_make_install ]];then

    echo -e " -- ${RED} ERROR: Failed to run make install for arpack-ng ${NC}"
    echo -e " -- ${RED} The error is: $fail_to_run_make_install ${NC}" 

    if [[ -d "$SUCCESS_FILE" ]]; then

	    rm $SUCCESS_FILE

    fi 
    
    exit 1 

  fi

  # Everything went fine, create the success file

  touch ../../arpack-ng-success.txt

fi
