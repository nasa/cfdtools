# This library is a copy of the xyzq_io library but compiled to support the
# extended versions of grid_block_structure defined in the tecplot_io library.
# In a perfect world, the tecplot version of these types would be extensions of
# the base type and the xzyq_io library would work with either version.
# However, the ability to do this sort generic programming is a recent addition
# to Fortran, so we would have to do an extensive rewrite of xyzq_io and client
# libraries to do it. Instead, we just compile mulitple versions of the
# library, since this is how CFDTOOLS has historicall solved the problem.

cfdtools_add_library(tecplot_xyzq_io
  # Use same sources as xyzq_io...
  ../xyzq_io/f_io.f90
  ../xyzq_io/xyq_io.f90
  ../xyzq_io/xyzq_io.f90
)
target_link_libraries(cfdtools_tecplot_xyzq_io
  PUBLIC
    # But link to tecplot-specific grid_block_structure
    cfdtools_tecplot_io
)
