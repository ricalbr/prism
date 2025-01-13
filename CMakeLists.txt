cmake_minimum_required(VERSION 3.12)
project(Prism)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-g -O3 -march=native -mavx2 -mfma -finline-limit=1000000 -ffp-contract=fast")

# Set output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Enable OpenMP for parallelism
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found.")
endif()

# Find YAML-CPP
find_package(yaml-cpp REQUIRED)

# Include directories
include_directories(include)

# Sources
file(GLOB_RECURSE SOURCES
    src/*.cpp
)

# Create the library
add_library(prism ${SOURCES})

# Link the library to YAML-CPP
target_link_libraries(prism PRIVATE yaml-cpp)

# Define the main executable
add_executable(prism_simulator main.cpp)

# Link the main executable to the library
target_link_libraries(prism_simulator PRIVATE prism)

# Installation rules
install(TARGETS prism_simulator DESTINATION bin)
install(DIRECTORY include/prism DESTINATION include)

