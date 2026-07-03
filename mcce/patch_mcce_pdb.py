#!/usr/bin/env python3
"""Patch a protein PDB using MCCE population and step2 conformer coordinates.

This script generates two PDBs:
- standard: patch titratable residues (and CYD mapped as CYS) with full atoms.
- amber: patch with amber-style residue-name conversions and limited hydrogens.
"""

from __future__ import annotations

import argparse
import bisect
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

TITRATABLE_BASES = {"ASP", "GLU", "LYS", "HIS", "CYS"}
POLAR_ROTAMER_BASES = {"ASN", "GLN"}
SPECIAL_BASES = {"CYD"}

# Backward-compatible aliases for callers that imported the old constants.
TARGET_BASES = TITRATABLE_BASES
TARGET_BASES_WITH_CYD = TITRATABLE_BASES | SPECIAL_BASES


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resid: int
    icode: str = ""


@dataclass
class AtomRecord:
    record: str
    serial: int
    atom_name: str
    resname: str
    chain: str
    resid: int
    icode: str
    x: float
    y: float
    z: float
    occ: float
    temp: float
    altloc: str = ""
    element: str = ""

    @property
    def residue_key(self) -> ResidueKey:
        return ResidueKey(self.chain, self.resid, self.icode)

    @property
    def identity_key(self) -> Tuple[str, int, str, str]:
        return (self.chain, self.resid, self.icode, self.atom_name.strip())


@dataclass(frozen=True)
class PopState:
    state_id: str
    populations: Tuple[float, ...]


@dataclass
class DominantChoice:
    base: str
    chain: str
    resid: int
    state_id: str
    population: float


POP_ID_RE = re.compile(r"^([A-Z0-9]{3})([+\-0-9A-Z]{2})(.)(\d{4})_(\d{3})$")
MCCE_RES_RE = re.compile(r"^(.)(\d{4})_(\d{3})$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Patch PDB by MCCE dominant states.")
    p.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    p.add_argument("--input-pdb", required=True, help="Original input PDB")
    p.add_argument("--population", required=True, help="MCCE population table")
    p.add_argument("--mcce-pdb", required=True, help="MCCE step2_out PDB")
    p.add_argument("--standard-out", default="patched_standard.pdb", help="Output standard PDB")
    p.add_argument("--amber-out", default="patched_amber.pdb", help="Output amber PDB")
    p.add_argument(
        "--report-out",
        default="population_report.txt",
        help="Population report text output",
    )
    p.add_argument(
        "--no-polar-rotamers",
        action="store_true",
        help="Do not patch ASN/GLN dominant rotamers",
    )
    return p.parse_args()


def infer_element(atom_name: str) -> str:
    core = atom_name.strip()
    if not core:
        return ""
    if core[0].isdigit():
        core = core[1:]
    if not core:
        return ""
    if len(core) >= 2 and core[1].islower():
        return core[:2].capitalize()
    return core[0].upper()


def format_atom_name(name: str) -> str:
    n = name.strip()
    if len(n) >= 4:
        return n[:4]
    if len(n) == 1:
        return f" {n}  "
    if len(n) == 2:
        if n[0].isdigit():
            return f"{n:<4}"
        return f" {n:<3}"
    if len(n) == 3:
        if n[0].isdigit():
            return f"{n:<4}"
        return f" {n}"
    return f"{n:>4}"


def parse_pdb_atom_line(line: str) -> AtomRecord:
    record = line[0:6].strip()
    return AtomRecord(
        record=record,
        serial=int(line[6:11]),
        atom_name=line[12:16].strip(),
        altloc=line[16:17].strip(),
        resname=line[17:20].strip(),
        chain=line[21:22],
        resid=int(line[22:26]),
        icode=line[26:27].strip(),
        x=float(line[30:38]),
        y=float(line[38:46]),
        z=float(line[46:54]),
        occ=float(line[54:60]) if line[54:60].strip() else 1.00,
        temp=float(line[60:66]) if line[60:66].strip() else 0.00,
        element=line[76:78].strip() if len(line) >= 78 else "",
    )


def atom_to_pdb_line(atom: AtomRecord, serial: int) -> str:
    atom_name = format_atom_name(atom.atom_name)
    element = atom.element.strip() or infer_element(atom.atom_name)
    return (
        f"{atom.record:<6}{serial:>5d} {atom_name}{(atom.altloc or ' '):1}"
        f"{atom.resname:>3} {(atom.chain or ' '):1}{atom.resid:>4d}{(atom.icode or ' '):1}   "
        f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occ:>6.2f}{atom.temp:>6.2f}"
        f"          {element:>2}"
    )


def ter_line(serial: int, resname: str, chain: str, resid: int, icode: str = "") -> str:
    return f"TER   {serial:>5d}      {resname:>3} {(chain or ' '):1}{resid:>4d}{(icode or ' '):1}"


def parse_original_pdb(path: Path):
    lines = path.read_text().splitlines()
    header: List[str] = []
    coord: List[Tuple[str, str]] = []
    conect: List[str] = []
    footer: List[str] = []

    in_coords = False
    in_conect = False

    for line in lines:
        rec = line[0:6].strip()
        if rec in {"ATOM", "HETATM", "TER"}:
            in_coords = True
            coord.append((rec, line))
            continue
        if rec == "CONECT":
            in_conect = True
            in_coords = False
            conect.append(line)
            continue
        if in_conect:
            footer.append(line)
        elif in_coords:
            coord.append((rec, line))
        else:
            header.append(line)
    return header, coord, conect, footer


def parse_population(
    path: Path,
    ph: float,
    target_bases: Iterable[str] = TARGET_BASES_WITH_CYD,
) -> Tuple[Dict[Tuple[str, str, int], DominantChoice], List[DominantChoice]]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError("Population file is empty")

    header_cols = lines[0].split()
    if not header_cols or header_cols[0].lower() != "ph":
        raise ValueError("Population header must start with 'ph'")

    ph_grid = [float(x) for x in header_cols[1:]]
    if ph < ph_grid[0] or ph > ph_grid[-1]:
        raise ValueError(f"Requested pH={ph} outside table range [{ph_grid[0]}, {ph_grid[-1]}]")

    target_base_set = set(target_bases)
    grouped: Dict[Tuple[str, str, int], List[PopState]] = {}

    for line in lines[1:]:
        cols = line.split()
        state_label = cols[0]
        m = POP_ID_RE.match(state_label)
        if not m:
            continue
        base, _two, chain, resid_s, state_id = m.groups()
        if base not in target_base_set:
            continue
        vals = tuple(float(x) for x in cols[1:])
        if len(vals) != len(ph_grid):
            raise ValueError(f"Population columns mismatch for {state_label}")
        key = (base, chain, int(resid_s))
        grouped.setdefault(key, []).append(PopState(state_id=state_id, populations=vals))

    def interpolate(vals: Sequence[float]) -> float:
        pos = bisect.bisect_left(ph_grid, ph)
        if pos < len(ph_grid) and ph_grid[pos] == ph:
            return vals[pos]
        i1 = pos
        i0 = pos - 1
        x0, x1 = ph_grid[i0], ph_grid[i1]
        y0, y1 = vals[i0], vals[i1]
        t = (ph - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)

    dominant_map: Dict[Tuple[str, str, int], DominantChoice] = {}
    dominant_list: List[DominantChoice] = []

    for (base, chain, resid), states in sorted(grouped.items(), key=lambda x: (x[0][1], x[0][2], x[0][0])):
        scored = [(st.state_id, interpolate(st.populations)) for st in states]
        scored.sort(key=lambda x: x[1], reverse=True)
        state_id, pop = scored[0]
        choice = DominantChoice(base=base, chain=chain, resid=resid, state_id=state_id, population=pop)
        dominant_map[(base, chain, resid)] = choice
        dominant_list.append(choice)

    return dominant_map, dominant_list


def parse_mcce_step2(path: Path) -> Dict[Tuple[str, str, int, str], List[AtomRecord]]:
    conf_atoms: Dict[Tuple[str, str, int, str], List[AtomRecord]] = {}

    with path.open() as f:
        for raw in f:
            if not raw.startswith(("ATOM", "HETATM")):
                continue
            # MCCE step2 lines are close to PDB fixed-width format. Parsing by
            # columns avoids failures when token and X coordinate touch, e.g.
            # "B0030_000-102.918".
            record = raw[0:6].strip()
            serial = int(raw[6:11])
            atom_name = raw[12:16].strip()
            resname = raw[17:20].strip()
            token = raw[21:30].strip()
            m = MCCE_RES_RE.match(token)
            if not m:
                continue
            chain, resid_s, state_id = m.groups()
            resid = int(resid_s)
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])

            atom = AtomRecord(
                record=record,
                serial=serial,
                atom_name=atom_name,
                resname=resname,
                chain=chain,
                resid=resid,
                icode="",
                x=x,
                y=y,
                z=z,
                occ=1.00,
                temp=0.00,
                altloc="",
                element=infer_element(atom_name),
            )
            conf_key = (resname, chain, resid, state_id)
            conf_atoms.setdefault(conf_key, []).append(atom)

    return conf_atoms


def select_atoms_for_residue(
    base: str,
    chain: str,
    resid: int,
    state_id: str,
    conf_atoms: Dict[Tuple[str, str, int, str], List[AtomRecord]],
) -> Tuple[List[AtomRecord], List[AtomRecord]]:
    bb_key = (base, chain, resid, "000")
    sc_key = (base, chain, resid, state_id)
    if bb_key not in conf_atoms:
        raise KeyError(f"Missing backbone conformer: {bb_key}")
    if sc_key not in conf_atoms:
        raise KeyError(f"Missing sidechain conformer: {sc_key}")
    return conf_atoms[bb_key], conf_atoms[sc_key]


def amber_resname(base: str, atom_names: Iterable[str]) -> str:
    names = set(atom_names)
    if base == "ASP":
        return "ASH" if ("HD1" in names or "HD2" in names) else "ASP"
    if base == "GLU":
        return "GLH" if ("HE1" in names or "HE2" in names) else "GLU"
    if base == "LYS":
        return "LYN" if ("HZ1" in names and "HZ2" in names and "HZ3" not in names) else "LYS"
    if base == "HIS":
        has_hd1 = "HD1" in names
        has_he2 = "HE2" in names
        if has_hd1 and has_he2:
            return "HIP"
        if has_hd1:
            return "HID"
        if has_he2:
            return "HIE"
        return "HIS"
    if base == "CYS":
        return "CYM" if "HG" not in names else "CYS"
    if base == "CYD":
        return "CYX"
    return base


def build_patch_definitions(
    dominant: Dict[Tuple[str, str, int], DominantChoice],
    conf_atoms: Dict[Tuple[str, str, int, str], List[AtomRecord]],
):
    std_patch: Dict[ResidueKey, List[AtomRecord]] = {}
    amb_patch: Dict[ResidueKey, List[AtomRecord]] = {}

    for (_, chain, resid), choice in sorted(dominant.items(), key=lambda x: (x[0][1], x[0][2], x[0][0])):
        base = choice.base
        bb_atoms, sc_atoms = select_atoms_for_residue(base, chain, resid, choice.state_id, conf_atoms)

        all_atoms = bb_atoms + sc_atoms

        # Standard: full atoms, CYD rendered as CYS in output.
        std_resname = "CYS" if base == "CYD" else base
        std_atoms = [
            AtomRecord(
                record="ATOM",
                serial=0,
                atom_name=a.atom_name,
                resname=std_resname,
                chain=chain,
                resid=resid,
                icode="",
                x=a.x,
                y=a.y,
                z=a.z,
                occ=1.00,
                temp=0.00,
                element=a.element,
            )
            for a in all_atoms
        ]
        std_patch[ResidueKey(chain, resid, "")] = std_atoms

        # Amber: keep only heavy atoms plus backbone H.
        amber_base_resname = amber_resname(base, (a.atom_name for a in sc_atoms))
        amber_atoms = []
        for a in all_atoms:
            is_heavy = not a.atom_name.startswith("H")
            keep = is_heavy or a.atom_name == "H"
            if not keep:
                continue
            amber_atoms.append(
                AtomRecord(
                    record="ATOM",
                    serial=0,
                    atom_name=a.atom_name,
                    resname=amber_base_resname,
                    chain=chain,
                    resid=resid,
                    icode="",
                    x=a.x,
                    y=a.y,
                    z=a.z,
                    occ=1.00,
                    temp=0.00,
                    element=a.element,
                )
            )
        amb_patch[ResidueKey(chain, resid, "")] = amber_atoms

    return std_patch, amb_patch


def patch_structure(
    coord_lines: Sequence[Tuple[str, str]],
    patch_map: Dict[ResidueKey, List[AtomRecord]],
) -> Tuple[List[Tuple[str, object]], Dict[int, Tuple[str, int, str, str]], Dict[Tuple[str, int, str, str], int]]:
    out_recs: List[Tuple[str, object]] = []

    # old serial -> identity key (for CONECT remap)
    old_id_map: Dict[int, Tuple[str, int, str, str]] = {}
    seen_patched: set[ResidueKey] = set()
    last_atom_for_ter: Optional[AtomRecord] = None

    for rec, line in coord_lines:
        if rec in {"ATOM", "HETATM"}:
            atom = parse_pdb_atom_line(line)
            old_id_map[atom.serial] = atom.identity_key
            rkey = atom.residue_key
            if rkey in patch_map:
                if rkey not in seen_patched:
                    for pa in patch_map[rkey]:
                        out_recs.append(("ATOM", pa))
                        last_atom_for_ter = pa
                    seen_patched.add(rkey)
                continue
            out_recs.append((rec, atom))
            last_atom_for_ter = atom
            continue

        if rec == "TER":
            if last_atom_for_ter is None:
                continue
            out_recs.append(("TER", last_atom_for_ter))
            continue

        out_recs.append(("RAW", line))

    # assign new serials and create output lines
    serial = 1
    new_lines: List[str] = []
    new_identity_map: Dict[Tuple[str, int, str, str], int] = {}

    for typ, obj in out_recs:
        if typ in {"ATOM", "HETATM"}:
            atom: AtomRecord = obj  # type: ignore[assignment]
            new_lines.append(atom_to_pdb_line(atom, serial))
            new_identity_map[atom.identity_key] = serial
            serial += 1
        elif typ == "TER":
            last_atom: AtomRecord = obj  # type: ignore[assignment]
            new_lines.append(ter_line(serial, last_atom.resname, last_atom.chain, last_atom.resid, last_atom.icode))
            serial += 1
        else:
            new_lines.append(str(obj))

    # Return the lines in tuple wrapper compatible with writer.
    wrapped = [("RAW", ln) for ln in new_lines]
    return wrapped, old_id_map, new_identity_map


def remap_conect_lines(
    conect_lines: Sequence[str],
    old_identity_by_serial: Dict[int, Tuple[str, int, str, str]],
    new_serial_by_identity: Dict[Tuple[str, int, str, str], int],
) -> List[str]:
    out: List[str] = []

    for line in conect_lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            src_old = int(parts[1])
            tgt_old = [int(x) for x in parts[2:]]
        except ValueError:
            continue

        src_id = old_identity_by_serial.get(src_old)
        if not src_id:
            continue
        src_new = new_serial_by_identity.get(src_id)
        if not src_new:
            continue

        mapped_targets: List[int] = []
        for old_t in tgt_old:
            t_id = old_identity_by_serial.get(old_t)
            if not t_id:
                continue
            t_new = new_serial_by_identity.get(t_id)
            if t_new and t_new not in mapped_targets:
                mapped_targets.append(t_new)

        if not mapped_targets:
            continue

        out.append("CONECT" + f"{src_new:>5d}" + "".join(f"{t:>5d}" for t in mapped_targets))

    return out


def write_pdb(path: Path, header: Sequence[str], coord_lines: Sequence[Tuple[str, str]], conect: Sequence[str], footer: Sequence[str]) -> None:
    out_lines: List[str] = []
    out_lines.extend(header)
    out_lines.extend([line for _, line in coord_lines])
    out_lines.extend(conect)
    out_lines.extend(footer)

    if not footer or (footer and footer[-1].strip() != "END"):
        out_lines.append("END")

    path.write_text("\n".join(out_lines) + "\n")


def format_report(
    ph: float,
    dominant_list: Sequence[DominantChoice],
    titratable_bases: Iterable[str] = TITRATABLE_BASES,
    polar_rotamer_bases: Iterable[str] = (),
) -> str:
    titratable_set = set(titratable_bases)
    polar_rotamer_set = set(polar_rotamer_bases)
    lines: List[str] = []
    lines.append(f"pH = {ph:.3f}")

    warn_count = 0

    def append_section(title: str, bases: set) -> None:
        nonlocal warn_count
        section_choices = [
            ch for ch in sorted(dominant_list, key=lambda x: (x.chain, x.resid, x.base)) if ch.base in bases
        ]
        if not section_choices:
            return
        if len(lines) > 1:
            lines.append("")
        lines.append(title)
        for ch in section_choices:
            pct = ch.population * 100.0
            warning = ""
            if 40.0 <= pct <= 60.0:
                warning = "  WARNING: dominant population in 40-60%"
                warn_count += 1
            lines.append(
                f"{ch.base} {ch.chain}{ch.resid:04d} state {ch.state_id}: {pct:6.2f}%{warning}"
            )

    append_section("Titratable residue dominant states (ASP/GLU/LYS/HIS/CYS):", titratable_set)
    append_section("Polar rotamer dominant states (ASN/GLN):", polar_rotamer_set)

    lines.append("")
    lines.append(f"Total warnings (40-60%): {warn_count}")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()

    input_pdb = Path(args.input_pdb)
    pop_path = Path(args.population)
    mcce_path = Path(args.mcce_pdb)
    standard_out = Path(args.standard_out)
    amber_out = Path(args.amber_out)
    report_out = Path(args.report_out)

    header, coord, conect, footer = parse_original_pdb(input_pdb)

    target_bases = TITRATABLE_BASES | SPECIAL_BASES
    polar_report_bases = set()
    if not args.no_polar_rotamers:
        target_bases |= POLAR_ROTAMER_BASES
        polar_report_bases = POLAR_ROTAMER_BASES

    dominant_map, dominant_list = parse_population(pop_path, args.ph, target_bases)
    conf_atoms = parse_mcce_step2(mcce_path)
    std_patch, amb_patch = build_patch_definitions(dominant_map, conf_atoms)

    std_coord, std_old_ids, std_new_ids = patch_structure(coord, std_patch)
    std_conect = remap_conect_lines(conect, std_old_ids, std_new_ids)
    write_pdb(standard_out, header, std_coord, std_conect, footer)

    amb_coord, amb_old_ids, amb_new_ids = patch_structure(coord, amb_patch)
    amb_conect = remap_conect_lines(conect, amb_old_ids, amb_new_ids)
    write_pdb(amber_out, header, amb_coord, amb_conect, footer)

    report_text = format_report(args.ph, dominant_list, TITRATABLE_BASES, polar_report_bases)
    report_out.write_text(report_text)
    print(report_text, end="")
    print(f"standard PDB written to: {standard_out}")
    print(f"amber PDB written to:    {amber_out}")
    print(f"population report:       {report_out}")


if __name__ == "__main__":
    main()
