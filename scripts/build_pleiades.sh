#!/bin/bash
# A script to build the "official" CFDTOOLS release for
# NASA's supercomputer, Pleiades. Run this script from
# the project root directory.
#
# Useage:
#   cd <dplr_root_directory>
#   ./scripts/build_pleiades.sh [toolchain [build_dir]]
#
# Arguments:
#   toolchain   Toolchain used to build CFDTOOLS: intel, gnu [def: intel]
#   build_dir   Directory where CFDTOOLS is built [def: build]
#               Final outputs will be in $build_dir/install

# Process arguments
toolchain="${1:-intel}"
build_dir_suffix=""
if [[ "${toolchain}" != "intel" ]]; then
    build_dir_suffix="_${toolchain}"
fi
build_dir="${2:-"build${build_dir_suffix}"}"

# Remember if script was sourced
sourced=""
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    sourced="true"
fi

# Configure the build environment
case $toolchain in
    intel)
        module purge
        module add pkgsrc/2021Q2
        module add comp-intel/2020.4.304
        export FC=ifort
        export CC=icc
        export CXX=icpc
        ;;
    gnu)
        module purge
        module add pkgsrc/2021Q2
        module add gcc/9.3
        export FC=gfortran
        export CC=gcc
        export CXX=g++
        ;;
    *)
        echo "ERROR: Unsupported toolchain (${toolchain})"
        if [[ ! ${sourced} ]]; then
            exit 1
        fi
        ;;
esac

# If we sourced the script, stop after configuring the env.
# If we executed the script, proceed with build
if [[ ${sourced} ]]; then
    echo "Environment configured for CFDTOOLS build ($toolchain)"
    return
fi

# Stop if any command below fails
set -e

# Create the build dir
if [ -d "$build_dir" ]; then
    echo "Previous build directory ($build_dir) exists! Delete? [y/n]"
    read response
    if [ "$response" == "y" ]; then
        rm -rf "$build_dir"
    else
        echo "Aborting"
        exit 1
    fi
fi
mkdir "$build_dir"
cd "$build_dir"

# Configure/build/install
cmake ..
make -j8 install
