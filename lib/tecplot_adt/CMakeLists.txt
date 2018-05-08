# This library is a copy of the adt library but compiled to support the
# extended versions of grid_block_structure defined in the tecplot_io library.
# In a perfect world, the tecplot version of these types would be extensions of
# the base type and the xzyq_io library would work with either version.
# However, the ability to do this sort generic programming is a recent addition
# to Fortran, so we would have to do an extensive rewrite of xyzq_io and client
# libraries to do it. Instead, we just compile mulitple versions of the
# library, since this is how CFDTOOLS has historicall solved the problem.

add_library(tecplot_adt
  ../adt/adt_utilities.f90
)
target_link_libraries(tecplot_adt
  PRIVATE
    global
  PUBLIC
    tecplot_io
)
target_include_directories(tecplot_adt PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
