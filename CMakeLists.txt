cmake_minimum_required(VERSION 3.0)

project(pygsw LANGUAGES Fortran)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_DEBUG_POSTFIX "" CACHE STRING "Debug file extension" FORCE)

add_subdirectory(extern/GSW-Fortran EXCLUDE_FROM_ALL)

add_library(pygsw_ SHARED pygsw/__init__.py pygsw/__main__.py src/wrapper.F90)
target_link_libraries(pygsw_ gsw_static)

# Find target Python interpreter.
find_package(PythonInterp)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pygsw" pygsw
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:pygsw_>" pygsw/
                   COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pygsw WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")