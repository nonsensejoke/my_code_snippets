"""
rechain_by_segi.py
------------------
Rewrite PDB chain IDs based on the segI field (columns 73-76).

Each unique segI value encountered (in order of first appearance)
is assigned a chain letter: A, B, C, D, ...

Usage:
    python rechain_by_segi.py input.pdb output.pdb

Only ATOM and HETATM records are modified (column 22 / chain ID).
All other lines are passed through unchanged.
"""

import sys
import string

def rechain(input_path: str, output_path: str) -> None:
    chain_pool = iter(string.ascii_uppercase)
    segi_to_chain: dict[str, str] = {}

    with open(input_path, "r") as f:
        lines = f.readlines()

    # First pass: build segI → chain letter mapping (order of first appearance)
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            segi = line[72:76].strip()
            if segi and segi not in segi_to_chain:
                segi_to_chain[segi] = next(chain_pool)

    print("segI → Chain ID mapping:")
    for segi, chain in segi_to_chain.items():
        print(f"  {segi:>4s}  →  {chain}")

    # Second pass: rewrite chain ID at column index 21 (1-based column 22)
    with open(output_path, "w") as out:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                segi = line[72:76].strip()
                chain = segi_to_chain.get(segi, line[21])
                line = line[:21] + chain + line[22:]
            out.write(line)

    print(f"\nOutput written to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rechain_by_segi.py input.pdb output.pdb")
        sys.exit(1)

    rechain(sys.argv[1], sys.argv[2])
