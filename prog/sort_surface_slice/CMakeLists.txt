add_executable(sort_surface_slice
  sort_surface_slice.f90
)
target_link_libraries(sort_surface_slice global)
target_link_libraries(sort_surface_slice
  numodules
  progtools
  ugridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME sort_surface_slice)
install(TARGETS sort_surface_slice DESTINATION ${CMAKE_INSTALL_BINDIR})
