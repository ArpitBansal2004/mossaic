# CMAKE generated file: DO NOT EDIT!
# Generated by CMake Version 3.22
cmake_policy(SET CMP0009 NEW)

# cs225_sources at lib/CMakeLists.txt:18 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/*.cpp")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/ColorSpace/ColorSpace.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/ColorSpace/Comparison.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/ColorSpace/Conversion.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/LUVAPixel.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/lib/cs225/PNG.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()

# lodepng_sources at lib/CMakeLists.txt:6 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/lib/lodepng/*.cpp")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/lib/lodepng/lodepng.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()

# lodepng_sources at lib/CMakeLists.txt:6 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/lib/lodepng/*.h")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/lib/lodepng/lodepng.h"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()

# util_sources at lib/CMakeLists.txt:12 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/lib/util/*.cpp")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/lib/util/coloredout.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/lib/util/util.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()

# src_sources at src/CMakeLists.txt:5 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/src/*.cpp")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/src/maptiles.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/src/mosaiccanvas.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/src/sourceimage.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/src/tileimage.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()

# tests_src at CMakeLists.txt:97 (file)
file(GLOB_RECURSE NEW_GLOB LIST_DIRECTORIES false "/workspaces/cs225/cs225git/mp_mosaics/tests/*.cpp")
set(OLD_GLOB
  "/workspaces/cs225/cs225git/mp_mosaics/tests/tests_part1.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/tests/tests_part2.cpp"
  "/workspaces/cs225/cs225git/mp_mosaics/tests/tests_student.cpp"
  )
if(NOT "${NEW_GLOB}" STREQUAL "${OLD_GLOB}")
  message("-- GLOB mismatch!")
  file(TOUCH_NOCREATE "/workspaces/cs225/cs225git/mp_mosaics/build/CMakeFiles/cmake.verify_globs")
endif()
