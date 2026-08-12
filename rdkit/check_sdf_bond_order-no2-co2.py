#!/usr/bin/env python
"""
Detect bond-order errors in SDF files whose first line holds the reference SMILES.

Verdict rule
------------
PASS/FAIL is decided by BOND ORDERS ONLY.

The discriminator is the per-atom sum of bond orders (explicit valence), compared
atom-by-atom against the reference SMILES:
  * it is invariant under the choice of Kekule structure, so a legitimately
    different but equivalent aromatic bond assignment is NOT flagged;
  * it is independent of formal charge, so files that merely omit the
    "M  CHG" block are NOT flagged.

Missing formal charges are reported as a NOTE, not an error, because that is a
file-writing convention issue rather than a wrong connectivity/bond order. To
prove it is harmless, the checker re-applies the expected charges and confirms
the mol block then reproduces the reference SMILES exactly.

Usage:  python check_bonds.py file1.sdf file2.sdf ...
Exit code is 1 if any file has a bond-order error.
"""

import sys
from collections import defaultdict

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BOND_ORDER = {
    Chem.BondType.SINGLE: 1.0,
    Chem.BondType.DOUBLE: 2.0,
    Chem.BondType.TRIPLE: 3.0,
    Chem.BondType.AROMATIC: 1.5,
}
NAME = {1.0: "single", 2.0: "double", 3.0: "triple", 1.5: "aromatic"}
MAXV = {"C": 4, "N": 3, "O": 2, "S": 6, "P": 5,
        "F": 1, "Cl": 1, "Br": 1, "I": 1, "B": 3}

GROUPS = [
    ("nitro -NO2", "[N+X3](=O)[O-]"),
    ("nitro -NO2", "[NX3](=O)=O"),
    ("carboxylate -COO-", "[CX3](=O)[O-]"),
    ("carboxyl -COOH", "[CX3](=O)[OX2H1]"),
]


def skeleton(mol):
    """Bond-order- and charge-stripped copy, for pure graph matching."""
    rw = Chem.RWMol(mol)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    return m


def valence_sums(mol):
    v = defaultdict(float)
    for b in mol.GetBonds():
        o = BOND_ORDER.get(b.GetBondType(), 0.0)
        v[b.GetBeginAtomIdx()] += o
        v[b.GetEndAtomIdx()] += o
    return v


def bond_orders(mol):
    return {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))):
            BOND_ORDER.get(b.GetBondType(), 0.0)
        for b in mol.GetBonds()
    }


def best_mapping(ref, sdf):
    """ref_idx -> sdf_idx mapping minimising bond-order conflicts."""
    matches = skeleton(sdf).GetSubstructMatches(
        skeleton(ref), uniquify=False, useChirality=False, maxMatches=20000
    )
    if not matches:
        return None
    rbo, sbo = bond_orders(ref), bond_orders(sdf)
    best, best_cost = None, None
    for m in matches:
        if any(ref.GetAtomWithIdx(ri).GetAtomicNum()
               != sdf.GetAtomWithIdx(si).GetAtomicNum()
               for ri, si in enumerate(m)):
            continue
        cost = sum(1 for (a, b), o in rbo.items()
                   if abs(sbo.get(tuple(sorted((m[a], m[b]))), -1) - o) > 1e-6)
        if best_cost is None or cost < best_cost:
            best, best_cost = m, cost
            if cost == 0:
                break
    return best


def group_labels(ref, r2s):
    labels = {}
    for label, sma in GROUPS:
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            continue
        for hit in ref.GetSubstructMatches(patt):
            for ri in hit:
                labels.setdefault(r2s[ri] + 1, label)
    return labels


def check(path):
    """Return True if the file passes the bond-order check."""
    name = path.split("/")[-1]
    text = open(path).read()
    smiles = text.splitlines()[0].strip()

    print("=" * 74)
    print(f"FILE   : {name}")
    print(f"SMILES : {smiles}")

    ref = Chem.MolFromSmiles(smiles)
    sdf = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
    if ref is None:
        print("VERDICT: ERROR - reference SMILES unparsable")
        return False
    if sdf is None:
        print("VERDICT: ERROR - mol block unparsable")
        return False
    ref = Chem.Mol(ref)
    Chem.Kekulize(ref, clearAromaticFlags=True)

    if ref.GetNumAtoms() != sdf.GetNumAtoms():
        print(f"VERDICT: ERROR - atom count {sdf.GetNumAtoms()} "
              f"vs {ref.GetNumAtoms()} in SMILES")
        return False

    m = best_mapping(ref, sdf)
    if m is None:
        print("VERDICT: ERROR - connectivity does not match the SMILES graph")
        return False
    r2s = {ri: si for ri, si in enumerate(m)}
    s2r = {si: ri for ri, si in r2s.items()}
    labels = group_labels(ref, r2s)

    vref, vsdf = valence_sums(ref), valence_sums(sdf)

    # ---- (1) bond orders: the only thing that decides the verdict ---------
    wrong_val = [(si, ri) for ri, si in r2s.items()
                 if abs(vsdf[si] - vref[ri]) > 1e-6]
    wrong_set = {si for si, _ in wrong_val}

    wrong_bonds, resonance = [], []
    for b in ref.GetBonds():
        ra, rb = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        sa, sb = r2s[ra], r2s[rb]
        sbond = sdf.GetBondBetweenAtoms(sa, sb)
        o_ref = BOND_ORDER.get(b.GetBondType(), 0.0)
        o_sdf = BOND_ORDER.get(sbond.GetBondType(), 0.0) if sbond else None
        if o_sdf is None or abs(o_sdf - o_ref) > 1e-6:
            rec = (min(sa, sb) + 1, max(sa, sb) + 1,
                   sdf.GetAtomWithIdx(sa).GetSymbol(),
                   sdf.GetAtomWithIdx(sb).GetSymbol(),
                   NAME.get(o_ref, o_ref), NAME.get(o_sdf, o_sdf))
            (wrong_bonds if (sa in wrong_set or sb in wrong_set)
             else resonance).append(rec)

    ok = not wrong_bonds and not wrong_val

    # ---- (2) formal charges: a note, never part of the verdict ------------
    missing_chg = [(si, ref.GetAtomWithIdx(ri).GetFormalCharge())
                   for ri, si in r2s.items()
                   if sdf.GetAtomWithIdx(si).GetFormalCharge()
                   != ref.GetAtomWithIdx(ri).GetFormalCharge()]

    print(f"VERDICT: {'OK - bond orders match the SMILES' if ok else 'BOND-ORDER ERROR'}")

    if wrong_bonds:
        print(f"\n  wrong bond orders ({len(wrong_bonds)}), SDF 1-based numbering:")
        for sa, sb, s1, s2, exp, got in sorted(wrong_bonds):
            grp = labels.get(sa) or labels.get(sb) or ""
            print(f"    bond {sa:>3}-{sb:<3} {s1}-{s2:<2} : SDF {got:<7}"
                  f" -> should be {exp:<7}" + (f"  [{grp}]" if grp else ""))

        print("\n  atoms with wrong valence sum (charge-independent):")
        for si, ri in sorted(wrong_val):
            grp = labels.get(si + 1, "")
            print(f"    atom {si+1:>3} {sdf.GetAtomWithIdx(si).GetSymbol():<2}"
                  f" valence {vsdf[si]:g}, expected {vref[ri]:g}"
                  + (f"  [{grp}]" if grp else ""))

        # over-valence, judged with the charges the SMILES calls for
        over = []
        for si in range(sdf.GetNumAtoms()):
            a = sdf.GetAtomWithIdx(si)
            q = ref.GetAtomWithIdx(s2r[si]).GetFormalCharge()
            if vsdf[si] > MAXV.get(a.GetSymbol(), 99) + q:
                over.append((si + 1, a.GetSymbol(), vsdf[si]))
        if over:
            print("\n  impossible valences (even after allowing the expected charge):")
            for i, sym, v in over:
                print(f"    atom {i:>3} {sym:<2} valence {v:g}")

    if resonance:
        print(f"\n  {len(resonance)} bond(s) differ only as Kekule/resonance "
              "bookkeeping - NOT an error:")
        for sa, sb, s1, s2, exp, got in sorted(resonance):
            print(f"    bond {sa:>3}-{sb:<3} {s1}-{s2:<2} : SDF {got}, SMILES {exp}")

    if missing_chg:
        print(f"\n  NOTE (not an error): 'M  CHG' block absent, "
              f"{len(missing_chg)} formal charge(s) dropped:")
        for si, q in sorted(missing_chg):
            print(f"    atom {si+1:>3} {sdf.GetAtomWithIdx(si).GetSymbol():<2}"
                  f" is 0 in the SDF, {q:+d} in the SMILES")
        if ok:
            # prove the charges are the ONLY thing missing
            rw = Chem.RWMol(sdf)
            for si, q in missing_chg:
                rw.GetAtomWithIdx(si).SetFormalCharge(q)
            fixed = rw.GetMol()
            try:
                Chem.SanitizeMol(fixed)
                same = (Chem.MolToSmiles(fixed)
                        == Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
            except Exception:
                same = False
            print("    -> restoring only these charges reproduces the reference "
                  "SMILES exactly: " + ("yes" if same else "no"))

    return ok


if __name__ == "__main__":
    results = [(p.split("/")[-1], check(p)) for p in sys.argv[1:]]
    if len(results) > 1:
        print("=" * 74)
        print("SUMMARY")
        for n, ok in results:
            print(f"  {'OK  ' if ok else 'FAIL'}  {n}")
    sys.exit(0 if all(ok for _, ok in results) else 1)
