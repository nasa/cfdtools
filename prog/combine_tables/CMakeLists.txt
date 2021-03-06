add_executable(combine_tables
  combine_tables.f
)
target_link_libraries(combine_tables cfdtools_global)
install(TARGETS combine_tables DESTINATION ${CMAKE_INSTALL_BINDIR})
