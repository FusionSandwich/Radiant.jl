#!/usr/bin/env bash
set -euo pipefail
g4env=/mnt/d/GEANT4/uw_hpge_geant4_model_final/.runtime/geant4-env
out=/mnt/c/rhts-proton/radiant/qualification/geant4_proton_compare
exe=/mnt/d/GEANT4/work/radiant_proton_compare_build/radiant_proton_compare
export PATH="$g4env/bin:$PATH"
export LD_LIBRARY_PATH="$g4env/lib:${LD_LIBRARY_PATH:-}"
export G4LEDATA="$g4env/share/Geant4/data/EMLOW8.8"
export G4ENSDFSTATEDATA="$g4env/share/Geant4/data/ENSDFSTATE3.0"
"$exe" "$out/geant4_stopping_low.csv" "$out/geant4_low_1p2_halfstep.csv" 200000 0.000005 1.2 0.00002
"$exe" "$out/geant4_stopping_low.csv" "$out/geant4_low_2_halfstep.csv" 200000 0.00001 2.0 0.00005
