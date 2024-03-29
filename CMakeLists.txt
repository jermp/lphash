cmake_minimum_required(VERSION 3.5)
project(LPHASH)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

MESSAGE(STATUS "Compiling for processor: " ${CMAKE_HOST_SYSTEM_PROCESSOR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  MESSAGE(STATUS "Compiling with flags: -march=native -mbmi2 -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
endif()

if (UNIX)

  MESSAGE(STATUS "Compiling with flags: -std=c++17 -O3 -ggdb -pthread -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # for multi-threading

  if (LPHASH_USE_SANITIZERS)
    MESSAGE(STATUS "Using sanitizers. Compiling with flags: -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()


set(LPHASH_SOURCES
  src/constants.cpp
  src/quartet_wtree.cpp
  src/minimizer.cpp
  src/partitioned_mphf.cpp
  src/unpartitioned_mphf.cpp
  src/mphf_utils.cpp
  src/parser_build.cpp
)

# main tool
add_executable(lphash src/lphash.cpp ${LPHASH_SOURCES})
target_link_libraries(lphash z)

# competitors
add_executable(ptbb_build src/ptbb/ptbb_build.cpp src/constants.cpp)
target_link_libraries(ptbb_build z)
add_executable(ptbb_query src/ptbb/ptbb_query.cpp src/constants.cpp)
target_link_libraries(ptbb_query z)

# tests
add_executable(test_sev tests/test_sev.cpp src/constants.cpp)
target_link_libraries(test_sev z)