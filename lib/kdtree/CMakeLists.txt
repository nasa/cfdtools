# Builds the fortran version only
add_library(kdtree
    src-f90/kdtree.f90
    src-f90/kdtree2.f90
)
target_link_libraries(kdtree PRIVATE cfdtools_global)
target_include_directories(kdtree PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES README DESTINATION doc RENAME kdtree)
