#!/usr/bin/env python3
"""
lipid21_to_martini3_popc.py
============================
Convert POPC from lipid21 atomistic PDB to Martini 3 coarse-grained .gro.

In lipid21, each POPC molecule is split into three residues:
    PA  - sn-1 palmitoyl chain  (C16:0)
    PC  - phosphocholine head + glycerol backbone
    OL  - sn-2 oleoyl chain     (C18:1, double bond C9=C10)

The script groups PA/PC/OL triplets, computes a center-of-geometry (COG)
for each Martini 3 bead, and writes a .gro file.

Mapping derived from:
    - CHARMM36 POPC topology (atom connectivity reference)
    - CHARMM -> Martini 2.2 map (vermouth/martinize2 source)
    - Martini 3 POPC.itp (bead names / types)

Key subtleties encoded in the mapping:
    - Glycerol C numbering is REVERSED between CHARMM and lipid21:
        lipid21 C3 (adj. phosphate)  = CHARMM C1  -> goes into GL1
        lipid21 C2 (sn-2 / oleoyl)  = CHARMM C2  -> goes into GL1
        lipid21 C1 (sn-1 / palmitoyl)= CHARMM C3 -> goes into GL2
    - GL2 spans two residues: C1 from PC + (O11, C11, O12) from PA
    - A chain (C1A..C4A) = sn-2 oleoyl (OL residue)
    - B chain (C1B..C4B) = sn-1 palmitoyl (PA residue)
    - Double bond C9=C10 (lipid21: C19=C110):
        C19  -> D2A  (last atom of that bead)
        C110 -> C3A  (first atom of that bead)

Requirements:
    pip install numpy

Usage:
    python lipid21_to_martini3_popc.py input.pdb output.gro
    python lipid21_to_martini3_popc.py input.pdb output.gro --box 15.0 15.0 10.0

PBC note:
    If your snapshot has molecules split across the periodic boundary,
    preprocess first:
        AMBER  -> cpptraj:  autoimage + trajout snapshot.pdb pdb
        GROMACS -> gmx trjconv -pbc whole -f traj.xtc -s topol.tpr
    The script applies a simple within-bead PBC unwrap when box dims
    are available, but whole-molecule unwrapping requires external tools.
"""

import sys
import argparse
import numpy as np
from collections import defaultdict


# ============================================================
#  lipid21 -> Martini 3 POPC atom-to-bead mapping
#  Format: { bead_name: { resname: [atom_names] } }
#  Only heavy atoms are listed (hydrogens auto-excluded).
# ============================================================

POPC_MAPPING = {
    "NC3": {"PC": ["N31",  "C31",  "C32",  "C33",  "C34",  "C35"]},
    "PO4": {"PC": ["P31",  "O31",  "O32",  "O33",  "O34"]},
    "GL1": {"PC": ["C3",   "C2",   "O21",  "C21",  "O22"]},
    "GL2": {"PC": ["C1"],
            "PA": ["O11",  "C11",  "O12"]},
    # --- sn-2 oleoyl chain (OL residue) ---
    "C1A": {"OL": ["C12",  "C13",  "C14",  "C15"]},
    "D2A": {"OL": ["C16",  "C17",  "C18",  "C19"]},           # C19 = double-bond C9
    "C3A": {"OL": ["C110", "C111", "C112", "C113", "C114"]},  # C110 = double-bond C10
    "C4A": {"OL": ["C115", "C116", "C117", "C118"]},
    # --- sn-1 palmitoyl chain (PA residue) ---
    "C1B": {"PA": ["C12",  "C13",  "C14",  "C15"]},
    "C2B": {"PA": ["C16",  "C17",  "C18",  "C19"]},
    "C3B": {"PA": ["C110", "C111", "C112", "C113"]},
    "C4B": {"PA": ["C114", "C115", "C116"]},
}

BEAD_ORDER = [
    "NC3", "PO4", "GL1", "GL2",
    "C1A", "D2A", "C3A", "C4A",
    "C1B", "C2B", "C3B", "C4B",
]


# ============================================================
#  PDB parser
# ============================================================

def parse_pdb(filename):
    """
    Read PDB file. Returns:
        residues  : { (resname, resnum): { atomname: np.array([x,y,z]) } }
        res_order : list of (resname, resnum) in order of first appearance
    Coordinates in Angstrom. Only PA/PC/OL residues are parsed.
    """
    residues  = {}
    res_order = []

    with open(filename) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue

            res_name  = line[17:20].strip()
            if res_name not in ("PA", "PC", "OL"):
                continue

            atom_name = line[12:16].strip()
            try:
                res_num = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            key = (res_name, res_num)
            if key not in residues:
                residues[key]  = {}
                res_order.append(key)

            residues[key][atom_name] = np.array([x, y, z])

    return residues, res_order


def get_box_from_cryst1(filename):
    """
    Extract orthorhombic box dimensions (nm) from CRYST1 record.
    Returns np.array([a,b,c]) or None.
    """
    with open(filename) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                try:
                    a = float(line[6:15]) / 10.0
                    b = float(line[15:24]) / 10.0
                    c = float(line[24:33]) / 10.0
                    return np.array([a, b, c])
                except ValueError:
                    pass
    return None


# ============================================================
#  Molecule grouping
# ============================================================

def group_popc_molecules(residues, res_order):
    """
    Group PA/PC/OL residues into POPC molecules.

    Assumption (standard lipid21 bilayer):
        - Each POPC = exactly one PA + one PC + one OL residue.
        - When the three residue types are each sorted by residue number,
          the i-th PA corresponds to the i-th PC and i-th OL.

    Returns list of dicts: [{"PA": atoms_dict, "PC": ..., "OL": ...}, ...]
    """
    def collect(rname):
        return sorted(
            [(rnum, residues[(rname, rnum)])
             for (rn, rnum) in res_order if rn == rname],
            key=lambda x: x[0]
        )

    pa_list = collect("PA")
    pc_list = collect("PC")
    ol_list = collect("OL")

    counts = {"PA": len(pa_list), "PC": len(pc_list), "OL": len(ol_list)}
    n = min(counts.values())

    if len(set(counts.values())) != 1:
        print(f"  Warning: residue counts differ {counts}. Using first {n} of each.")

    return [
        {"PA": pa_list[i][1],
         "PC": pc_list[i][1],
         "OL": ol_list[i][1]}
        for i in range(n)
    ]


# ============================================================
#  CG bead position
# ============================================================

def unwrap_positions(positions, box_ang):
    """
    Apply minimum-image convention within a set of positions
    (used to fix within-bead PBC artifacts).
    Uses positions[0] as reference. Modifies in place.

    positions : list of np.array([x,y,z]) in Angstrom
    box_ang   : np.array([bx,by,bz]) in Angstrom, or None to skip
    """
    if box_ang is None or len(positions) <= 1:
        return positions

    ref = positions[0].copy()
    result = [ref]
    for pos in positions[1:]:
        delta = pos - ref
        delta -= np.round(delta / box_ang) * box_ang
        result.append(ref + delta)
    return result


def compute_bead_cog(mol, bead_name, box_ang=None):
    """
    Center of geometry (COG) for one CG bead.

    mol       : { "PA": {atomname: pos}, "PC": ..., "OL": ... }
    bead_name : Martini bead name
    box_ang   : box dimensions in Angstrom for within-bead PBC fix (or None)

    Returns np.array([x,y,z]) in Angstrom, or None if no atoms found.
    """
    positions = []
    missing   = []

    for res_name, atom_names in POPC_MAPPING[bead_name].items():
        res_atoms = mol.get(res_name, {})
        for aname in atom_names:
            if aname in res_atoms:
                positions.append(res_atoms[aname])
            else:
                missing.append(f"{res_name}/{aname}")

    if missing:
        print(f"    [!] {bead_name}: atoms not found: {', '.join(missing)}")

    if not positions:
        return None

    positions = unwrap_positions(positions, box_ang)
    return np.mean(positions, axis=0)


# ============================================================
#  .gro writer
# ============================================================

def write_gro(filename, cg_molecules, box_nm):
    """
    Write GROMACS .gro file.

    cg_molecules : list of { bead_name -> np.array([x,y,z]) in Angstrom }
    box_nm       : np.array([bx,by,bz]) in nm
    """
    total_atoms = len(cg_molecules) * len(BEAD_ORDER)

    with open(filename, "w") as fh:
        fh.write("POPC Martini 3 CG - converted from lipid21\n")
        fh.write(f"{total_atoms:5d}\n")

        atom_num = 1
        for mol_idx, bead_dict in enumerate(cg_molecules):
            res_num = (mol_idx % 99999) + 1   # .gro max 5 digits

            for bead_name in BEAD_ORDER:
                pos = bead_dict.get(bead_name)
                if pos is None:
                    pos = np.zeros(3)

                x_nm, y_nm, z_nm = pos / 10.0  # Å → nm

                fh.write(
                    f"{res_num:5d}"
                    f"{'POPC':<5s}"
                    f"{bead_name:>5s}"
                    f"{(atom_num % 100000):5d}"
                    f"{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
                )
                atom_num += 1

        # Box vector line (nm, orthorhombic)
        fh.write(f"   {box_nm[0]:.5f}   {box_nm[1]:.5f}   {box_nm[2]:.5f}\n")


# ============================================================
#  Validation helper
# ============================================================

def validate_molecule(mol_idx, bead_dict):
    """
    Sanity-check one CG molecule.
    Returns True if all beads present and head-tail distance is plausible.
    """
    ok = True

    # All beads must exist
    for bead in BEAD_ORDER:
        if bead_dict.get(bead) is None:
            print(f"  [!] Mol {mol_idx+1}: bead {bead} is missing.")
            ok = False

    if ok:
        # NC3 (head) - C4B (tail) distance should be ~2-5 nm in a bilayer
        head = bead_dict["NC3"]
        tail = bead_dict["C4B"]
        dist_nm = np.linalg.norm(head - tail) / 10.0
        if not (0.5 < dist_nm < 6.0):
            print(f"  [!] Mol {mol_idx+1}: NC3-C4B distance = {dist_nm:.2f} nm "
                  f"(expected 1-4 nm). PBC artifact?")
            ok = False

    return ok


# ============================================================
#  Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Convert lipid21 POPC → Martini 3 CG .gro",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Requirements:")[0].strip()
    )
    parser.add_argument("input",  help="Input PDB file (lipid21 POPC)")
    parser.add_argument("output", help="Output .gro file (Martini 3 CG)")
    parser.add_argument(
        "--box", nargs=3, type=float, metavar=("X", "Y", "Z"),
        help="Box dimensions in nm (overrides CRYST1 record)"
    )
    parser.add_argument(
        "--no-pbc-fix", action="store_true",
        help="Skip within-bead PBC unwrapping"
    )
    args = parser.parse_args()

    # ── Parse input ───────────────────────────────────────────
    print(f"Reading:  {args.input}")
    residues, res_order = parse_pdb(args.input)

    pa_n = sum(1 for rn, _ in res_order if rn == "PA")
    pc_n = sum(1 for rn, _ in res_order if rn == "PC")
    ol_n = sum(1 for rn, _ in res_order if rn == "OL")
    print(f"Residues: {pa_n} PA  {pc_n} PC  {ol_n} OL")

    # ── Box dimensions ────────────────────────────────────────
    if args.box:
        box_nm = np.array(args.box)
        print(f"Box (arg): {box_nm[0]:.3f} x {box_nm[1]:.3f} x {box_nm[2]:.3f} nm")
    else:
        box_nm = get_box_from_cryst1(args.input)
        if box_nm is not None:
            print(f"Box (CRYST1): {box_nm[0]:.3f} x {box_nm[1]:.3f} x {box_nm[2]:.3f} nm")
        else:
            print("Warning: no box dimensions found.")
            print("         Provide with --box X Y Z, or add CRYST1 to PDB.")
            print("         Using placeholder 10 x 10 x 10 nm.")
            box_nm = np.array([10.0, 10.0, 10.0])

    box_ang = box_nm * 10.0 if not args.no_pbc_fix else None

    # ── Group molecules ───────────────────────────────────────
    print("Grouping PA/PC/OL triplets into POPC molecules ...")
    molecules = group_popc_molecules(residues, res_order)
    print(f"POPC molecules: {len(molecules)}")

    # ── Compute CG beads ──────────────────────────────────────
    print("Computing COG for each bead ...")
    cg_molecules = []
    n_warn = 0

    for i, mol in enumerate(molecules):
        bead_dict = {
            bead: compute_bead_cog(mol, bead, box_ang)
            for bead in BEAD_ORDER
        }
        if not validate_molecule(i, bead_dict):
            n_warn += 1
        cg_molecules.append(bead_dict)

    # ── Write output ──────────────────────────────────────────
    print(f"Writing:  {args.output}")
    write_gro(args.output, cg_molecules, box_nm)

    total_beads = len(molecules) * len(BEAD_ORDER)
    print()
    print(f"Done:  {len(molecules)} POPC  x  12 beads  =  {total_beads} CG atoms")
    if n_warn:
        print(f"       {n_warn} molecules had warnings - inspect carefully")

    print()
    print("Next steps:")
    print("  1. Merge CG POPC with other system components")
    print("     (protein martinize2, water/ions from insane or CHARMM-GUI)")
    print("  2. Use Martini 3 POPC.itp for bonded parameters")
    print("  3. Energy minimization  →  NVT  →  NPT  →  production MD")
    print()
    print("If you see PBC artifacts (NC3-C4B dist > 4 nm), preprocess with:")
    print("  AMBER   : cpptraj  → autoimage → trajout snapshot.pdb pdb")
    print("  GROMACS : gmx trjconv -pbc whole -f in.xtc -s topol.tpr -o out.pdb")


if __name__ == "__main__":
    main()