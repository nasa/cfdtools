cfdtools_add_executable(Stardust_Lines Stardust_Lines.f90)
target_link_libraries(cfdtools_Stardust_Lines
  cfdtools_adt
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_interp3d
  cfdtools_surface_patch_utilities
  cfdtools_xyzq_io
)

cfdtools_add_executable(Stardust_Integration Stardust_Integration.f90)
target_link_libraries(cfdtools_Stardust_Integration
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_integratelib
  cfdtools_interp3d
  cfdtools_intrinsics
  cfdtools_searchlib
  cfdtools_xyzq_io
)

cfdtools_add_executable(point_source point_source.f90)




