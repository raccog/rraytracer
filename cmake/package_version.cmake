cmake_minimum_required(VERSION 3.7...3.21)

# Use git to retrieve version
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    OUTPUT_VARIABLE RRAYTRACER_VERSION
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
    get_target_property(RRAYTRACER_VERSION rraytracer VERSION)
endif()

# Split version
string(REPLACE "." ";" RRAYTRACER_VERSION_LIST ${RRAYTRACER_VERSION})
list(GET RRAYTRACER_VERSION_LIST 0 RRAYTRACER_VERSION_MAJOR)
list(GET RRAYTRACER_VERSION_LIST 1 RRAYTRACER_VERSION_MINOR)
list(GET RRAYTRACER_VERSION_LIST 2 RRAYTRACER_VERSION_PATCH)
