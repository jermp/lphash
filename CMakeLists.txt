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

if (UNIX AND NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
endif()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lz")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")        # for multi-threading
  if (LPHASH_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
endif()

# set(Z_LIB_SOURCES
#   include/gz/zip_stream.cpp
# )

set(LPHASH_SOURCES
  src/constants.cpp
  src/parser_build.cpp
  src/parser_query.cpp
  src/quartet_wtree.cpp
  src/mphf.cpp
)

# Create a static lib
# add_library(lphash_static STATIC
#   ${Z_LIB_SOURCES}
#   ${LPHASH_SOURCES}
# )

# build with the lib
add_executable(build src/build.cpp ${LPHASH_SOURCES})
add_executable(query src/query.cpp ${LPHASH_SOURCES})

# target_link_libraries(build
#   lphash_static
#   z
# )

# add_executable(bench src/bench.cpp ${SSHASH_SOURCES})

# add_executable(check src/check.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
# target_link_libraries(check z)

# add_executable(query src/query.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
# target_link_libraries(query z)

# add_executable(permute src/permute.cpp ${SSHASH_SOURCES} ${Z_LIB_SOURCES})
# target_link_libraries(permute z)