add_executable(sixseries
  phep.f
  sixseries.f
)
target_link_libraries(sixseries cfdtools_global)
install(TARGETS sixseries DESTINATION ${CMAKE_INSTALL_BINDIR})
