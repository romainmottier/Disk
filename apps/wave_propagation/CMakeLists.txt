set(LINK_LIBS diskpp)

# Multithreaded Pardiso: diskpp links MKL with the Intel threading layer (mkl_intel_thread + libiomp5). Where libiomp5
# is LLVM's libomp (Ubuntu packages), the Pardiso factorization crashes unless MKL_NUM_THREADS=1. With GCC, the GNU
# threading layer and libgomp are linked first so that MKL uses them (mkl_gnu_thread is taken next to mkl_core, another
# MKL version may be installed elsewhere).
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND MKL_mkl_core_LIBRARY AND OpenMP_gomp_LIBRARY)
    get_filename_component(MKL_LIBRARY_DIR ${MKL_mkl_core_LIBRARY} DIRECTORY)
    find_library(MKL_GNU_THREAD_LIBRARY NAMES mkl_gnu_thread PATHS ${MKL_LIBRARY_DIR} NO_DEFAULT_PATH)
    if (MKL_GNU_THREAD_LIBRARY)
        set(LINK_LIBS ${MKL_GNU_THREAD_LIBRARY} ${OpenMP_gomp_LIBRARY} ${LINK_LIBS})
    endif()
endif()

# The two prototypes (first- and second-order formulations), configured by a Lua input file
add_executable(wave_propagation src/wave_propagation.cpp)
target_link_libraries(wave_propagation ${LINK_LIBS})

# Former driver with all the research prototypes (select one in its main()), not built by default
option(WAVE_PROPAGATION_LEGACY "Build wave_propagation_legacy, the former driver with all the research prototypes" OFF)
if (WAVE_PROPAGATION_LEGACY)
    add_executable(wave_propagation_legacy src/wave_propagation_legacy.cpp)
    target_link_libraries(wave_propagation_legacy ${LINK_LIBS})
endif()

# Lua input files (annotated templates and examples) and the reproducibility campaigns of earlier papers
# (share/reproducibility/<papers>/<paper>: python3 run_campaign.py in each paper folder) copied next to the
# executable at configure time: edit the copies under another name, or the sources in share/. The campaigns write
# their results and figures in the build copies only, kept by later configurations.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/share/first_order.lua ${CMAKE_CURRENT_SOURCE_DIR}/share/second_order.lua
          ${CMAKE_CURRENT_SOURCE_DIR}/share/examples ${CMAKE_CURRENT_SOURCE_DIR}/share/reproducibility
          ${CMAKE_CURRENT_SOURCE_DIR}/share/PARAMETRES.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# location of diskpp's meshes for the campaign (polygonal meshes of the FVCA5 benchmark)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/reproducibility/papers_Duran/paths.json
     "{ \"meshes\": \"${PROJECT_SOURCE_DIR}/refactor_old_diskpp_code/meshes\" }\n")
# marker of the build copy of the papers_Mottier campaigns (results are written there, not in the sources)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/reproducibility/papers_Mottier/paths.json
     "{ \"sources\": \"${CMAKE_CURRENT_SOURCE_DIR}/share/reproducibility/papers_Mottier\" }\n")

install(TARGETS wave_propagation RUNTIME DESTINATION bin)
