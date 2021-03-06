# This library is a copy of the surface_patch_utilities library but compiled to
# support the extended versions of grid_block_structure defined in the
# tecplot_io library.  In a perfect world, the tecplot version of these types
# would be extensions of the base type and this library would work with either
# version.  However, the ability to do this sort generic programming is a
# recent addition to Fortran, so we would have to do an extensive rewrite of
# the library and the client programs to do it. Instead, we just compile
# mulitple versions of the library, since this is how CFDTOOLS has historicall
# solved the problem.

add_library(tecplot_surface_patch_utilities
  ../surface_patch_utilities/surface_patch_utilities.f90
  ../surface_patch_utilities/triangulate_patches.f90
)
target_link_libraries(tecplot_surface_patch_utilities
  PRIVATE
    cfdtools_global
  PUBLIC
    geomlib
    linsys
    tecplot_io
    triangulation_io
    numodules
)
target_include_directories(tecplot_surface_patch_utilities PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
