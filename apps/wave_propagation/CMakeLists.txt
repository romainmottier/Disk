execute_process(
    COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(pybind11_DIR ${pybind11_dir})
find_package(pybind11 REQUIRED)

########################################################

set(LINK_LIBS diskpp)

add_executable(wave_propagation src/wave_propagation.cpp)
target_link_libraries(wave_propagation ${LINK_LIBS} pybind11::embed)

install(TARGETS wave_propagation RUNTIME DESTINATION bin)
