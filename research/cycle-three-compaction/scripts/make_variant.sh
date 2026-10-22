#!/bin/bash
# Build the compact cycle-three variant from the accepted proof plus the files in ../lean.
# usage: scripts/make_variant.sh <repo-root> <output-folder>
set -euo pipefail
root=$1; out=$2; here=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$out"
cp "$root"/data/cycle-three/{*.lean,proof.json,solution.md} "$out"/
# New syzygy and differential certificates (drop-in replacements).
cp "$here"/lean/NineKronecker.lean "$here"/lean/FullSyzygy.lean "$here"/lean/CubicDifferential.lean "$out"/
# The only consumer of the old syzygy layer: one import and one theorem name.
sed -i 's/^import Proposal.CubicSyzygy$/import Proposal.FullSyzygy/; s/CycleCovariants.IdentityData.full_syzygy/CycleCovariants.Full.full_syzygy/' "$out"/PhysicalCubicJets.lean
# Modules no longer reachable from Solution.lean (CubicDominance was already unused).
( cd "$out" && rm -f BivariateExpressions.lean CubicCoefficients{A,B,C,D,E,F,G,H}.lean CubicDominance.lean \
    CubicEncodingCheck.lean CubicGradeCheck.lean CubicIdentityData.lean CubicIdentityEvaluation.lean \
    CubicIdentityProducts.lean CubicKronecker.lean CubicSyzygy.lean )
echo "variant written to $out: $(ls "$out"/*.lean | wc -l) Lean files, $(cat "$out"/*.lean | wc -c) bytes"
