# Version range of CMake required
cmake_minimum_required(VERSION 3.7...3.21)

# Load rtga library
add_subdirectory(rtga)

# Include header files of dependencies
target_include_directories(rraytracer PUBLIC rtga/include)
