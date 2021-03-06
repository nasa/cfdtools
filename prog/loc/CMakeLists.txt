add_executable(loc
  loc.f
)
target_link_libraries(loc cfdtools_global)
target_link_libraries(loc
  progtools
)
install(FILES README DESTINATION doc RENAME loc)
install(TARGETS loc DESTINATION ${CMAKE_INSTALL_BINDIR})
