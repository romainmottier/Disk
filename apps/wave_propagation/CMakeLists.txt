set(LINK_LIBS diskpp)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

add_executable(wave_propagation src/wave_propagation.cpp)

target_link_libraries(wave_propagation
    PRIVATE
    ${LINK_LIBS}
    pybind11::embed
    Python::Python
)

install(TARGETS wave_propagation RUNTIME DESTINATION bin)

