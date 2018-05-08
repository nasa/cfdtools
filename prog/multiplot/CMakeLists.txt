add_executable(multiplot
  bounds.f
  intrinsics.f
  multiplot.f
  psutil.f
)
target_link_libraries(multiplot global intrinsics ugridlib)
install(TARGETS multiplot DESTINATION ${CMAKE_INSTALL_BINDIR})
