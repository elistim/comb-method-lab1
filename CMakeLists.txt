cmake_minimum_required(VERSION 3.10)
project(my_program C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Set compiler warnings for the active toolchain
if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall)
endif()

# Set output directory for built executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.c
    src/app.c
    src/tasks/task_6_1/task_6_1.c
    src/tasks/task_6_1/expression_6_1.c
    src/tasks/task_6_1/function_provider_6_1.c
    src/tasks/task_6_1/solver_6_1.c
    src/tasks/task_6_2/task_6_2_1.c
    src/tasks/task_6_2/task_6_2_2.c
    src/tasks/task_6_2/task_6_2_3.c
)

# Create executable
add_executable(my_program ${SOURCES})
