#!/usr/bin/env bash
set -euo pipefail

PYMOL="${PYMOL:-pymol}"

for cif in "$@"; do
    [ -f "$cif" ] || { echo "skip (not found): $cif" >&2; continue; }
    cif_abs="$(realpath "$cif")"
    pdb_abs="${cif_abs%.cif}.pdb"
    pml="${cif_abs%.cif}.pml"

    cat > "$pml" <<EOF
load $cif_abs
set retain_order, 1
save $pdb_abs
EOF

    "$PYMOL" -cqk "$pml"
    rm -f "$pml"
    echo "✓ $cif -> $pdb_abs"
done
