find_package(OpenMP)

if(OpenMP_CXX_FOUND AND USVG_USE_OPENMP)
  target_link_libraries(${LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()