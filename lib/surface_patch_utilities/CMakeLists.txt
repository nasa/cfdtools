add_library(surface_patch_utilities
  surface_grid_structure.f90.unused
  surface_neighbors.f90.incomplete
  surface_patch_utilities.f90
  surface_patch_utilities.f90.before-surface-tangent
  triangulate_patches.f90
)
target_link_libraries(surface_patch_utilities global)
install(FILES README DESTINATION doc)
