# Defines any CFDTOOLS-specific compilations flags needed in addition to
# the default flags provided by CMake. For reference, the default CMake
# flags (same for C/C++/Fortran) are:
#
#   Debug:          -g
#   Release:        -O3 -DNDEBUG
#   RelWithDebInfo: -O2 -DNDEBUG -g
#   MinSizeRel:     -Os -DNDEBUG

set(options "")

#-----------------------------------------------------------------------
# Intel
#-----------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    # Common
    list(APPEND options
        -real-size 64       # Default to 64bit reals
        -fpp                # Run preprocessor on all files
        -extend-source 132  # Allows long lines in source files
        -assume buffered_io # FCONVERT is hopelessly slow otherwise
    )

    # Debug
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        list(APPEND options
            -check all      # Runtime bounds checking
            -fpe0           # Abort if get NaN, etc.
            -traceback      # Emit stacktrace on error
            -warn all       # Full warning analysis
        )
    endif()


#-----------------------------------------------------------------------
# GNU
#-----------------------------------------------------------------------
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

    # Common
    list(APPEND options
        -cpp                        # Run preprocessor on all files
        -fdefault-real-8            # Default to 64bit reals
        -fdefault-double-8          # Default to 64bit doubles
        -ffixed-line-length-none    # Allow long lines in source files
        -ffree-line-length-none
    )

    # Debug
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        list(APPEND options
            -fcheck=all                       # Runtime bounds checking
            -ffpe-trap=invalid,zero,overflow  # Abort if NaN, etc.
            -fbacktrace                       # Emit stacktrace on error
            -Wall -Wextra                     # Full warning analysis
        )
    endif()


#-----------------------------------------------------------------------
# Unknown Compiler
#-----------------------------------------------------------------------
else()
    message(
        FATAL_ERROR
        "Fortran compiler '${CMAKE_Fortran_COMPILER_ID}' is not currently supported."
    )
endif()


#-----------------------------------------------------------------------
# Cache Results
#-----------------------------------------------------------------------
set(CFDTOOLS_COMPILE_OPTIONS ${options} CACHE STRING
    "Compiler options added to all DPLR components"
    FORCE)
mark_as_advanced(CFDTOOLS_COMPILE_OPTIONS)
