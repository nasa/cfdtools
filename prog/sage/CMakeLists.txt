add_executable(sage
  numutil.f
  sage.f
)
target_link_libraries(sage cfdtools_global)
install(TARGETS sage DESTINATION ${CMAKE_INSTALL_BINDIR})
