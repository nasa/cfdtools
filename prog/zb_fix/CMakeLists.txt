add_executable(zb_fix
  zb-fix.f90
)
target_link_libraries(zb_fix cfdtools_global)
install(TARGETS zb_fix DESTINATION ${CMAKE_INSTALL_BINDIR})
