cmake_minimum_required(VERSION 3.15)

project(pygsw LANGUAGES Fortran C)

if (MSVC)
  add_compile_options(/libs:static)
endif()

set(CMAKE_DEBUG_POSTFIX "" CACHE STRING "Debug file extension" FORCE)

add_subdirectory(extern/GSW-Fortran EXCLUDE_FROM_ALL)

add_library(cgsw pygsw/__init__.py pygsw/__main__.py src/wrapper.F90)
target_link_libraries(cgsw gsw_static)

# Find target Python interpreter.
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT _pygsw.c
  COMMENT "Creating ${CMAKE_CURRENT_BINARY_DIR}/_pygsw.c from ${CMAKE_CURRENT_SOURCE_DIR}/src/_pygsw.pyx"
  COMMAND Python3::Interpreter -m cython "${CMAKE_CURRENT_SOURCE_DIR}/src/_pygsw.pyx" -o "${CMAKE_CURRENT_BINARY_DIR}/_pygsw.c"
  DEPENDS src/_pygsw.pyx
  VERBATIM
)

Python3_add_library(_pygsw MODULE ${CMAKE_CURRENT_BINARY_DIR}/_pygsw.c)
target_link_libraries(_pygsw PRIVATE cgsw Python3::NumPy)

# Build wheel
add_custom_command(OUTPUT pygsw_wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pygsw" pygsw
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:_pygsw>" pygsw/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(pygsw_wheel ALL DEPENDS pygsw_wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --force-reinstall --upgrade pygsw WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")
