cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD_REQUIRED 20)

# Set Build Type
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    message("Building Release")
else()
    message("Building Debug")
endif()

# Set -A Platform
if(WIN32)
    set(CMAKE_SYSTEM_PROCESSOR win32)
    message("Configuration set: win32")
else()
    message("Default configuration")
endif()

# Set -T Toolset
# if(WIN32)
#     set(CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE x64)
#     message("Toolset set: x64")
# else()
#     message("Default configuration")
# endif()

#==============================

set(PROJECT_NAME Analyze)
project(${PROJECT_NAME})

find_package(ROOT REQUIRED)
set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

add_subdirectory(src/vendor/catima)
include_directories(src)

file(GLOB SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} catima)

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR})