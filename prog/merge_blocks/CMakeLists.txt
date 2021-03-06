add_executable(merge_blocks
  merge_blocks.f90
)
add_executable(merge_blocks_with_halos
  merge_blocks_with_halos.f90
)
target_link_libraries(merge_blocks cfdtools_global progtools)
target_link_libraries(merge_blocks_with_halos cfdtools_global progtools)
install(FILES README DESTINATION doc RENAME merge_blocks)
install(TARGETS merge_blocks merge_blocks_with_halos DESTINATION ${CMAKE_INSTALL_BINDIR})
