add_executable(merge_tables
  merge_tables.f90
)
target_link_libraries(merge_tables cfdtools_global)
target_link_libraries(merge_tables
  progtools
)
install(FILES README DESTINATION doc RENAME merge_tables)
install(TARGETS merge_tables DESTINATION ${CMAKE_INSTALL_BINDIR})
