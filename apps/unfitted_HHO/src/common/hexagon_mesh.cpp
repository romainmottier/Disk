/* Hexagonal meshes of the FVCA5 benchmark (the mesher of diskpp defines
 * non-inline functions: it is included in this translation unit only) */
#include "diskpp/mesh/meshgen.hpp"
#include "meshes.hpp"

namespace unfitted_meshes {

void make_hexagon_mesh(disk::generic_mesh<double,2>& msh, size_t level)
{
    auto mesher = disk::make_fvca5_hex_mesher(msh);
    mesher.make_level(level);
}

} // namespace unfitted_meshes
