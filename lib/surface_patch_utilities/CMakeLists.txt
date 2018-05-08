add_library(surface_patch_utilities
  surface_patch_utilities.f90
  triangulate_patches.f90
)
target_link_libraries(surface_patch_utilities
  PRIVATE
    global
  PUBLIC
    grid_block_structure
)
target_include_directories(surface_patch_utilities PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME surface_patch_utilities)
