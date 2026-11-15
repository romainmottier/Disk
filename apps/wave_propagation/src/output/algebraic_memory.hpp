/*
 *       /\        Wave propagation (DiskPP)
 *
 * Algebraic memory of a run (output.algebraic_memory, opt-in): the memory of the linear algebra of the time
 * integration only, without the libraries, the mesh, the assembly and the outputs:
 *   operator  the evolution operator used by the time steps (matrix-free blocks, or the row-major sparse blocks of the
 *             fast mode of the second-order formulation), counted from its containers;
 *   state     the solution vector;
 *   scheme    what the time scheme allocates and keeps during the integration: its work vectors, the levels of MLTS
 *             (coefficients, chains, lists of unknowns, local operators), the fine system of LI and its factors. It is
 *             the largest growth of the heap in use (glibc mallinfo2: blocks allocated and not freed, whatever the
 *             allocator keeps or returns to the system) at the end of the steps (probe), above its value just before
 *             the integration (baseline). The forcing is evaluated once before the baseline, so that its caches are
 *             not counted.
 * Deterministic: it does not depend on the resident set size, on the threads or on the duration of the run.
 */

#ifndef wave_propagation_algebraic_memory_hpp
#define wave_propagation_algebraic_memory_hpp

#include <malloc.h>
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace algebraic_memory {

    struct state_type {
        bool on = false;
        size_t base = 0, peak = 0;     // heap in use before the integration, largest value at the probes
        size_t op = 0, x = 0;          // operator and state, counted from their containers
    };

    inline state_type& S() {
        static state_type s;
        return s;
    }

    // bytes allocated by malloc and not freed (all the arenas, and the blocks obtained by mmap)
    inline size_t heap_in_use() {
        const struct mallinfo2 m = mallinfo2();
        return m.uordblks + m.hblkhd;
    }

    inline void enable(bool on) { S().on = on; }
    inline bool enabled() { return S().on; }

    // just before the time integration: operator and state bytes, and the heap in use
    inline void baseline(size_t operator_bytes, size_t state_bytes) {
        if (!S().on) return;
        S().op = operator_bytes;
        S().x = state_bytes;
        S().base = S().peak = heap_in_use();
    }

    // at the end of a time step (no-op unless enabled)
    inline void probe() {
        if (!S().on) return;
        S().peak = std::max(S().peak, heap_in_use());
    }

    inline size_t scheme_bytes() { return S().peak > S().base ? S().peak - S().base : 0; }

    // "algebraic memory: operator X MB, state Y MB, time scheme Z MB, total T MB" (MB = 2^20 bytes)
    inline std::string report() {
        const double mb = 1024.0 * 1024.0;
        const size_t z = scheme_bytes();
        std::ostringstream os;
        os << std::fixed << std::setprecision(3) << "      algebraic memory: operator " << double(S().op) / mb << " MB, state "
           << double(S().x) / mb << " MB, time scheme " << double(z) / mb << " MB, total " << double(S().op + S().x + z) / mb
           << " MB";
        return os.str();
    }

}

#endif /* wave_propagation_algebraic_memory_hpp */
