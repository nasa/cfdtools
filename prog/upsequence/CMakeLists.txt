add_executable(upsequence
  upsequence.f90
)
target_link_libraries(upsequence global)
target_link_libraries(upsequence
  xyzq_io
)
install(FILES README DESTINATION doc RENAME upsequence)
install(TARGETS upsequence DESTINATION ${CMAKE_INSTALL_BINDIR})
