add_executable(extract_peaks
  extract_peaks.f90
)
target_link_libraries(extract_peaks global)
install(FILES README DESTINATION doc RENAME extract_peaks)
install(TARGETS extract_peaks DESTINATION ${CMAKE_INSTALL_BINDIR})
