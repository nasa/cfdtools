#!/bin/bash
# A script to build the "official" CFDTOOLS release for
# NASA's supercomputer, Pleiades. Run this script from
# the project root directory.
#
# Useage:
#   cd <dplr_root_directory>
#   ./scripts/build_pleiades.sh [preset]
#
# Arguments:
#   preset      CMake preset used for build: intel, gnu, etc. [def: intel]

# Process arguments
preset="${1:-intel}"

# Remember if script was sourced
sourced=""
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    sourced="true"
fi

# Configure the build environment
case "${preset}" in
    *intel|deploy-nas)
        module purge
        module add pkgsrc/2021Q2
        module add comp-intel/2020.4.304
        ;;
    *gnu)
        module purge
        module add pkgsrc/2021Q2
        module add gcc/9.3
        ;;
    *)
        echo "ERROR: No module configuration specifed for preset ${preset}"
        if [[ ! ${sourced} ]]; then
            exit 1
        fi
        ;;
esac

# If we sourced the script, stop after configuring the env.
# If we executed the script, proceed with build
if [[ ${sourced} ]]; then
    echo "Environment configured for CFDTOOLS development (${preset})"
    return
fi

# Configure/build/install
cmake --preset="${preset}" -G Ninja
cmake --build "build-${preset}" --target install -j8

