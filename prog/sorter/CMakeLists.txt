add_executable(sorter
  sorter.f
)
target_link_libraries(sorter cfdtools_global searchlib progtools)
install(TARGETS sorter DESTINATION ${CMAKE_INSTALL_BINDIR})
