add_executable(compare_blocks
  compare_blocks.f
)
target_link_libraries(compare_blocks global)
install(FILES README DESTINATION doc RENAME compare_blocks)
install(TARGETS compare_blocks DESTINATION ${CMAKE_INSTALL_BINDIR})
