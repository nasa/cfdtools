add_library(adt
  adt_utilities.f90
)
target_link_libraries(adt
  PRIVATE
    cfdtools_global
    geomlib
  PUBLIC
    grid_block_structure
)
target_include_directories(adt PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME adt)
