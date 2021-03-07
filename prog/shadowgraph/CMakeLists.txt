# This library uses a custom version of grid_block_structure (grrr...). As
# such, we cannot link against the standard xyzq_io, and grid_block_utilities
# libraries by specifying targets like we normally would. So instead we
# explicitly recompile the sources that we need from those libraries. In a
# perfect world, the custom grid_block_type would be extensions of basic type
# and the xzyq_io and grid_block_utilities libraries would work with either
# version.  However, the ability to do this sort generic programming is a
# recent addition to Fortran, so we would have to do an extensive rewrite of
# xyzq_io and client libraries to do it. Did I mention I love fortran?
cfdtools_add_executable(shadowgraph
  FieldView_Data.f90
  flow_gradients.f90
  grid_block_structure.f90
  process_hex_cell.f90
  process_prism.f90
  process_pyramid.f90
  process_tet_cell.f90
  read_fieldview.f90
  shadowgraph.f90
  triangulate_quad_face.f90
  VTK_Data.f90
  # Recompile libraries for custom grid_block_structure
  ../../lib/xyzq_io/xyzq_io.f90
  ../../lib/grid_block_utilities/grid_block_utilities.f90
)
target_link_libraries(cfdtools_shadowgraph
  cfdtools_interp3d
  cfdtools_intrinsics
  cfdtools_numodules
  cfdtools_progtools
)
