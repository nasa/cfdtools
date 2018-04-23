add_executable(compare_patches
  compare_patches.f
)
target_link_libraries(compare_patches global)
install(FILES README DESTINATION doc RENAME compare_patches)
install(TARGETS compare_patches DESTINATION ${CMAKE_INSTALL_BINDIR})
