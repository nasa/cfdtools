cfdtools_add_executable(gestimate
  gestimate.f90
)
target_link_libraries(cfdtools_gestimate cfdtools_linsys)
