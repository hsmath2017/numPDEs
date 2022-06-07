add_subdirectory(Core)
add_subdirectory(MultiGrid)

add_library(${MAIN_LIB}
    ${Core_Files}
    ${MultiGrid_Files}
)
# target_include_directories(${MAIN_LIB} PUBLIC ./ ${CMAKE_SOURCE_DIR}/example/common)
#target_link_libraries(${MAIN_LIB} PUBLIC hdf5-shared siloh5)
# target_include_directories(${MAIN_LIB} PUBLIC ${HDF5_INCLUDE_DIRS} ${SILO_INCLUDE_DIRS})
# target_link_directories(${MAIN_LIB} PUBLIC ${HDF5_LIBRARY_DIRS} ${SILO_LIBRARY_DIRS})
# target_link_libraries(${MAIN_LIB} PUBLIC ${HDF5_LIBRARIES} ${SILO_LIBRARIES})

if(USE_MKL)
  if(APPLE)
    target_link_libraries(${MAIN_LIB} PUBLIC -L${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
  elseif(WIN32)
    message(SEND_ERROR "Using Intel(R) MKL on Windows system is not supported yet. ")
  else()
    target_link_libraries(${MAIN_LIB} PUBLIC -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl)
  endif()
elseif(USE_AOCL)
  target_include_directories(${MAIN_LIB} PUBLIC ${AOCL_ROOT}/include)
  target_link_directories(${MAIN_LIB} PUBLIC ${AOCL_ROOT}/lib)
  target_link_libraries(${MAIN_LIB} PUBLIC lapacke flame blis-mt -no-pie)
elseif(USE_OPENBLAS)
  target_include_directories(${MAIN_LIB} PUBLIC ${OpenBLAS_INCLUDE_DIRS})
  target_link_libraries(${MAIN_LIB} PUBLIC ${OpenBLAS_LIBRARIES} pthread gfortran)
else()
  target_include_directories(${MAIN_LIB} PUBLIC ${BLAS_INCLUDE_DIRS} ${CBLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${LAPACKE_INCLUDE_DIRS})
  target_link_directories(${MAIN_LIB} PUBLIC ${BLAS_LIBRARY_DIRS} ${CBLAS_LIBRARY_DIRS} ${LAPACK_LIBRARY_DIRS} ${LAPACKE_LIBRARY_DIRS})
  target_link_libraries(${MAIN_LIB} PUBLIC ${BLAS_LIBRARIES} ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES} gfortran)
endif()


