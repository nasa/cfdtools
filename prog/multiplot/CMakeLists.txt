add_executable(multiplot
  bounds.f
  intrinsics.f
  multiplot.f
  psutil.f
)
target_link_libraries(multiplot cfdtools_global intrinsics ugridlib)
install(TARGETS multiplot DESTINATION ${CMAKE_INSTALL_BINDIR})
