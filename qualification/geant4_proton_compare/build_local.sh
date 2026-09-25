#!/usr/bin/env bash
set -euo pipefail

g4env=/mnt/d/GEANT4/uw_hpge_geant4_model_final/.runtime/geant4-env
source_dir=/mnt/c/rhts-proton/radiant/qualification/geant4_proton_compare
build_dir=/mnt/d/GEANT4/work/radiant_proton_compare_build
export PATH="$g4env/bin:$PATH"
export LD_LIBRARY_PATH="$g4env/lib:${LD_LIBRARY_PATH:-}"
export G4LEDATA="$g4env/share/Geant4/data/EMLOW8.8"

"$g4env/bin/cmake" -S "$source_dir" -B "$build_dir" -G Ninja \
  -DGeant4_DIR="$g4env/lib/cmake/Geant4" \
  -DCMAKE_PREFIX_PATH="$g4env" \
  -DEXPAT_INCLUDE_DIR="$g4env/include" \
  -DEXPAT_LIBRARY="$g4env/lib/libexpat.so" \
  -DCMAKE_CXX_COMPILER="$g4env/bin/x86_64-conda-linux-gnu-c++" \
  -DCMAKE_MAKE_PROGRAM="$g4env/bin/ninja" \
  -DCMAKE_BUILD_TYPE=Release
"$g4env/bin/cmake" --build "$build_dir" --parallel 1
