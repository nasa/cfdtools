add_executable(merge_files
  merge_files.f90
)
target_link_libraries(merge_files cfdtools_global)
target_link_libraries(merge_files
  tecplot_io
  progtools
)
install(FILES README DESTINATION doc RENAME merge_files)
install(TARGETS merge_files DESTINATION ${CMAKE_INSTALL_BINDIR})
