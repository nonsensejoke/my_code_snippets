"""
SMARTS for detecting secondary amide (peptide) bonds with RDKit.

Use case: identify amide/peptide bonds so their central N-C bond can be
excluded from the list of freely rotatable bonds (peptide bonds are planar
and rigid due to resonance, so they should NOT be treated as free dihedrals).

    amide_smarts = '[H]-[N;H1]-[CX3]=[OX1]'

Atom-by-atom meaning:
    [H]      explicit hydrogen on the amide nitrogen (requires removeHs=False)
    [N;H1]   nitrogen bearing exactly one hydrogen (secondary amide N)
    [CX3]    carbonyl carbon, exactly 3 connections
    [OX1]    carbonyl oxygen, single connection (the C=O oxygen)

Why the leading [H] matters
---------------------------
The leading [H] is a REAL atom node, so it occupies index 0 of every match
tuple. The resulting match ordering is:

    index 0 -> H     (amide hydrogen)
    index 1 -> N     (nitrogen)
    index 2 -> C     (carbonyl carbon)
    index 3 -> O     (carbonyl oxygen)

That is convenient: the central peptide-bond atoms are match[1] (N) and
match[2] (C). If you drop the leading [H], indices shift left by one and the
N-C pair becomes match[0], match[1].

Note: [H] is redundant for *selecting* which nitrogens match, because [N;H1]
already requires exactly one H. Its real job here is index alignment / making
the H explicit in the match.

Gotchas / limitations
---------------------
- Requires explicit hydrogens: read the molecule with removeHs=False, otherwise
  the [H] node never matches.
- Silently MISSES tertiary / N-substituted amides because [N;H1] demands exactly
  one H. The most important case is the PROLINE peptide bond (its backbone N has
  no H) - such bonds will NOT be flagged and may be treated as rotatable.
- Also misses primary amides C(=O)NH2 ([N;H2]) - usually fine, since those are
  not backbone peptide bonds.
- This is a *secondary amide* detector, not a strict peptide-bond detector: any
  non-peptide secondary amide in the molecule will also match (alpha carbons are
  not constrained).

RDKit substructure matches do not guarantee bond-atom ordering, so when removing
the matched N-C pair from a rotatable-bond list, check BOTH orderings:

    if (b[1], b[2]) in rotatable: rotatable.remove((b[1], b[2]))
    if (b[2], b[1]) in rotatable: rotatable.remove((b[2], b[1]))

A broader alternative that also covers proline / tertiary amides:

    amide_smarts_any = '[NX3]-[CX3]=[OX1]'   # N-C match at index 0, 1
"""

from rdkit import Chem


def find_amide_bonds(mol):
    """Return substructure matches for secondary amide bonds.

    mol must be built with explicit hydrogens (removeHs=False).
    Each match is a 4-tuple (H_idx, N_idx, C_idx, O_idx).
    """
    amide_smarts = '[H]-[N;H1]-[CX3]=[OX1]'
    patt = Chem.MolFromSmarts(amide_smarts)
    return mol.GetSubstructMatches(patt)


def strip_amide_bonds_from_rotatable(rotatable_bonds, amide_matches):
    """Remove the central N-C bond of each amide from a rotatable-bond list.

    rotatable_bonds: iterable of (a, b) atom-index tuples
    amide_matches:   output of find_amide_bonds (H, N, C, O per match)
    """
    rot = list(rotatable_bonds)
    for _h, n, c, _o in amide_matches:
        if (n, c) in rot:
            rot.remove((n, c))
        if (c, n) in rot:
            rot.remove((c, n))
    return rot
