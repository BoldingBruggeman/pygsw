cmake_minimum_required(VERSION 3.0)

project(pygsw LANGUAGES Fortran)

set(CMAKE_DEBUG_POSTFIX "" CACHE STRING "Debug file extension" FORCE)

add_subdirectory(extern/GSW-Fortran EXCLUDE_FROM_ALL)

add_library(cgsw pygsw/__init__.py pygsw/__main__.py src/wrapper.F90)
target_link_libraries(cgsw gsw_static)

# Find target Python interpreter.
set(Python_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/_pygsw.c
  COMMAND ${Python3_EXECUTABLE} -m cython src/_pygsw.pyx
  MAIN_DEPENDENCY src/_pygsw.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

Python3_add_library(_pygsw MODULE src/_pygsw.c)
target_link_libraries(_pygsw PRIVATE cgsw Python3::NumPy)

# Build wheel
add_custom_command(OUTPUT pygsw_wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pygsw" pygsw
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_pygsw>" pygsw/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(pygsw_wheel ALL DEPENDS pygsw_wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pygsw WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
