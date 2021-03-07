add_executable(sorter
  sorter.f
)
target_link_libraries(sorter cfdtools_global cfdtools_searchlib cfdtools_progtools)
install(TARGETS sorter DESTINATION ${CMAKE_INSTALL_BINDIR})
