set( THE_TARGET "radix")
message(STATUS "Configuration of ${THE_TARGET}")
project(${THE_TARGET})
cmake_minimum_required(VERSION 2.8)

set(LIBRARY_OUTPUT_PATH     "${CMAKE_BINARY_DIR}/lib" CACHE PATH "Output directory for libraries")
set(EXECUTABLE_OUTPUT_PATH  "${CMAKE_BINARY_DIR}/bin" CACHE PATH "Output directory for applications")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

set(TNNS_CMAKE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${TNNS_CMAKE_DIR}")

#find_package(TBB REQUIRED)
#include_directories(${TBB_INCLUDE_DIRS})

#find_package(OpenMP REQUIRED)

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

#set(CMAKE_C_FLAGS "-std=gnu++0x")
#set(CMAKE_CXX_FLAGS "-std=gnu++0x")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT")

#include("${TNNS_CMAKE_DIR}/DetectIntelCompiler.cmake")
#find_intel_compiler()

#set(INTEL_SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/with_cilk.cpp")
#set(CILK_RADIX_OBJ "${CMAKE_BINARY_DIR}/with_cilk.obj")

#add_custom_command(OUTPUT ${CILK_RADIX_OBJ}
#                   COMMAND ${INTEL_COMPILER_EXE}
#                   ARGS /c /Iinclude /I${INTEL_INCLUDE_DIR} /Fo${CILK_RADIX_OBJ} ${INTEL_SOURCE_FILE}
#                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#                   MAIN_DEPENDENCY ${INTEL_SOURCE_FILE})

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_20,code=compute_20")
cuda_compile(CUDA_RADIX_OBJ  "${CMAKE_CURRENT_SOURCE_DIR}/src/radix.cu")
message("!!!!!!11${CUDA_RADIX_OBJ}")

#set(CMAKE_C_FLAGS "-std=gnu++11")
#set(CMAKE_CXX_FLAGS "-std=gnu++11")

add_executable( ${THE_TARGET} "${CMAKE_CURRENT_SOURCE_DIR}/src/radix.cpp"
                              "${CMAKE_CURRENT_SOURCE_DIR}/src/sequential.cpp"
                              "${CMAKE_CURRENT_SOURCE_DIR}/src/timing.cpp"
                              #"${CMAKE_CURRENT_SOURCE_DIR}/src/with_tbb.cpp"
                              #"${CMAKE_CURRENT_SOURCE_DIR}/src/with_omp.cpp"
                              "${CMAKE_CURRENT_SOURCE_DIR}/src/with_cuda.cpp"
                               ${CUDA_RADIX_OBJ}
                              #"${CILK_RADIX_OBJ}"
                            )

set(${THE_TARGET}_LINK_LIST "")

#if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
#    list(APPEND ${THE_TARGET}_LINK_LIST ${TBB_DEBUG_LIBRARIES})
#    set(TBB_SO_NAME "tbb_debug")
#else()
#    list(APPEND ${THE_TARGET}_LINK_LIST ${TBB_LIBRARIES})
#    set(TBB_SO_NAME "tbb")
#endif()

#list(APPEND ${THE_TARGET}_LINK_LIST "${INTEL_LIBRARY_PATH}/libmmt.lib"
#                                    "${INTEL_LIBRARY_PATH}/libirc.lib"
#                                    "${INTEL_LIBRARY_PATH}/svml_dispmt.lib"
#                                    "${INTEL_LIBRARY_PATH}/libdecimal.lib"
#                                    "${INTEL_LIBRARY_PATH}/cilkrts.lib")

target_link_libraries(${THE_TARGET} ${${THE_TARGET}_LINK_LIST} ${CUDA_LIBRARIES})

set_target_properties(${THE_TARGET} PROPERTIES
        DEBUG_POSTFIX "D"
        ARCHIVE_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
        LINKER_LANGUAGE CXX
    )

#string(REPLACE "${TBB_INSTALL_DIR}/lib/" "/" TBB_PLATFORM_COMPILER_PATH ${TBB_LIBRARY_DIRS})
#set(TBB_SO_PATH "${TBB_INSTALL_DIR}/bin${TBB_PLATFORM_COMPILER_PATH}")
#add_custom_command(TARGET ${THE_TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND}
#-E copy "${TBB_SO_PATH}/${TBB_SO_NAME}.dll" "${EXECUTABLE_OUTPUT_PATH}/${TBB_SO_NAME}.dll")
