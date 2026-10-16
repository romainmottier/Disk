#!/bin/bash
# two-size meshes of E2 (in ./meshes): 16 x 16 squares, the central 2 x 2 split down to the depth L (abrupt transition)
G=$(cd "$(dirname "$0")/../../../../src/mesh_generation/lts_efficiency" && pwd)/generate_refined_quadtree.py
mkdir -p meshes
for L in 1 2 3 4 5; do
  python3 $G --n 16 --depth $L --grading two --feature box:0.4375,0.5625,0.4375,0.5625 --output meshes/two_n16_m2_L$L.txt
done
