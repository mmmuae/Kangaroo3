cmake_minimum_required(VERSION 3.20)
project(Kangaroo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(APPLE)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Architectures to build" FORCE)

  # Maximum Apple Silicon M1/M2/M3 optimizations
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=apple-m1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native")

  # Aggressive optimization flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvectorize")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fslp-vectorize")

  # Link-Time Optimization
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

  # Cache line alignment for M1 Max (64 bytes)
  add_compile_definitions(CACHE_LINE_SIZE=64)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG" CACHE STRING "" FORCE)

set(KANGAROO_SOURCES
  Backup.cpp
  Check.cpp
  HashTable.cpp
  Kangaroo.cpp
  Merge.cpp
  Network.cpp
  PartMerge.cpp
  Thread.cpp
  Timer.cpp
  main.cpp
  SECPK1/Int.cpp
  SECPK1/IntGroup.cpp
  SECPK1/IntMod.cpp
  SECPK1/Point.cpp
  SECPK1/Random.cpp
  SECPK1/SECP256K1.cpp
)

add_executable(kangaroo ${KANGAROO_SOURCES})
target_include_directories(kangaroo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(kangaroo PRIVATE -Wall -Wextra -Wpedantic)
endif()

find_package(Threads REQUIRED)
target_link_libraries(kangaroo PRIVATE Threads::Threads)

install(TARGETS kangaroo RUNTIME DESTINATION bin)
