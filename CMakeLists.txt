#------------------------------------------------------------------------------
# Preamble
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.3)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CFDTOOLS_VERSION_MAJOR  15)
set(CFDTOOLS_VERSION_MINOR  0)
set(CFDTOOLS_VERSION_PATCH  0)
set(CFDTOOLS_VERSION "${CFDTOOLS_VERSION_MAJOR}.${CFDTOOLS_VERSION_MINOR}.${CFDTOOLS_VERSION_PATCH}")

project(cfdtools
    VERSION ${CFDTOOLS_VERSION}
    LANGUAGES Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


#------------------------------------------------------------------------------
# Global Defaults
#------------------------------------------------------------------------------
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install"
        CACHE PATH "DPLR installation directory"
        FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' since none specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(compile_options) # After set CMAKE_BUILD_TYPE
add_library(global INTERFACE)
target_compile_options(global INTERFACE ${CFDTOOLS_COMPILE_OPTIONS})


#------------------------------------------------------------------------------
# Include Project Subdirectories
#------------------------------------------------------------------------------
add_subdirectory(lib/aa_nurbs)
add_subdirectory(lib/adt)
add_subdirectory(lib/blas2lib)
add_subdirectory(lib/blaslib)
add_subdirectory(lib/dt_nurbs)
add_subdirectory(lib/eigenlib)
add_subdirectory(lib/forsythe)
add_subdirectory(lib/geomlib)
add_subdirectory(lib/grid_block_structure)
add_subdirectory(lib/grid_block_utilities)
add_subdirectory(lib/gridlib)
add_subdirectory(lib/grmodules)
add_subdirectory(lib/integrate)
add_subdirectory(lib/interp2d)
add_subdirectory(lib/interp3d)
add_subdirectory(lib/interplib)
add_subdirectory(lib/kdtree)
add_subdirectory(lib/linsys)
add_subdirectory(lib/numodules)
add_subdirectory(lib/obsolete)
add_subdirectory(lib/optlib)
add_subdirectory(lib/progtools)
add_subdirectory(lib/searchlib)
add_subdirectory(lib/special_functions)
add_subdirectory(lib/surface_patch_utilities)
add_subdirectory(lib/table_io)
add_subdirectory(lib/tecplot_io)
add_subdirectory(lib/triangulation_io)
add_subdirectory(lib/ugridlib)
add_subdirectory(lib/xyzq_io)

