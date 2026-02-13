cmake_minimum_required(VERSION 3.16)
project(wave_propagation CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development)

execute_process(
    COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT PYBIND11_DIR)
    message(FATAL_ERROR "pybind11 not found! Please install it in your Python environment.")
endif()

list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")

find_package(pybind11 REQUIRED)

set(LINK_LIBS diskpp)

add_executable(wave_propagation src/wave_propagation.cpp)

target_link_libraries(wave_propagation
    PRIVATE
    ${LINK_LIBS}
    pybind11::embed
    Python::Python
)

install(TARGETS wave_propagation RUNTIME DESTINATION bin)

