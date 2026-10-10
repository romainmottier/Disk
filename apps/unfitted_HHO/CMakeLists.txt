set(LINK_LIBS diskpp)

# one translation unit per type of mesh (compiled in parallel)
add_executable(unfitted_HHO src/unfitted_HHO.cpp
               src/run_quadrangles.cpp src/run_triangles.cpp src/run_polygons.cpp
               src/run_hexahedra.cpp src/run_tetrahedra.cpp src/run_polyhedra.cpp
               src/common/hexagon_mesh.cpp)

# The Intel threading layer of MKL linked by diskpp (mkl_intel_thread) needs
# the OpenMP runtime of Intel. Where libiomp5 is the libomp of LLVM (Debian,
# Ubuntu), the multithreaded factorization of Pardiso fails. With GCC, the GNU
# threading layer of MKL and libgomp are linked first: their symbols take
# precedence over those of the Intel layer (the layer is taken from the
# directory of the MKL core library found by diskpp).
if (HAVE_INTEL_MKL AND MKL_mkl_core_LIBRARY AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_filename_component(MKL_CORE_DIR ${MKL_mkl_core_LIBRARY} DIRECTORY)
    find_library(MKL_GNU_THREAD_LIBRARY NAMES mkl_gnu_thread HINTS ${MKL_CORE_DIR} NO_DEFAULT_PATH)
    find_package(OpenMP COMPONENTS CXX)
    if (MKL_GNU_THREAD_LIBRARY AND OpenMP_gomp_LIBRARY)
        message(STATUS "unfitted_HHO: GNU threading layer of MKL")
        set(LINK_LIBS ${MKL_GNU_THREAD_LIBRARY} ${OpenMP_gomp_LIBRARY} ${LINK_LIBS})
    endif()
endif()

target_link_libraries(unfitted_HHO ${LINK_LIBS})
install(TARGETS unfitted_HHO RUNTIME DESTINATION bin)
