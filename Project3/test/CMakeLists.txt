# configure_file("${CMAKE_SOURCE_DIR}/test/dirConfig.h.in" 
# "${CMAKE_SOURCE_DIR}/test/dirConfig.h")
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
# Tests on the rectangular domains
add_executable(testVCycle testVCycle.cpp)
target_link_libraries(testVCycle ${MAIN_LIB})
add_executable(testOp testOperator.cpp)
target_link_libraries(testOp ${MAIN_LIB})