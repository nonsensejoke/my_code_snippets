#!/usr/bin/env python3
"""Export all MCCE microstates for one residue into a PyMOL session."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


MCCE_RES_RE = re.compile(r"^(.)(\d{4})_(\d{3})$")
FORT38_LABEL_RE = re.compile(r"^([A-Z0-9]{3})([+\-0-9A-Z]{2})(.)(\d{4})_(\d{3})$")
TARGET_RE = re.compile(r"^([A-Za-z])(\d{1,4})$")
AXIS_TO_VECTOR = {
    "x": (1.0, 0.0, 0.0),
    "y": (0.0, 1.0, 0.0),
    "z": (0.0, 0.0, 1.0),
}


@dataclass(frozen=True)
class AtomRecord:
    record: str
    atom_name: str
    resname: str
    chain: str
    resid: int
    state_id: str
    conf_code: str
    x: float
    y: float
    z: float
    element: str


@dataclass(frozen=True)
class Step2Data:
    atoms_by_conformer: Dict[Tuple[str, str, int, str], List[AtomRecord]]
    conf_codes: Dict[Tuple[str, str, int, str], str]


@dataclass(frozen=True)
class Microstate:
    label: str
    object_name: str
    resname: str
    chain: str
    resid: int
    state_id: str
    atoms: Tuple[AtomRecord, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a PyMOL .pse showing all MCCE microstates for one residue."
    )
    parser.add_argument("--step2-pdb", required=True, help="MCCE step2_out.pdb file")
    parser.add_argument("--residue", required=True, help="Target residue ID, e.g. B0071")
    parser.add_argument(
        "--fort38",
        default=None,
        help="MCCE fort.38 population table; defaults to fort.38 in the current directory, if present",
    )
    parser.add_argument("--out-prefix", default=None, help="Output prefix for .pse/.pml/temp PDBs")
    parser.add_argument("--spacing", type=float, default=6.0, help="Spacing between microstates in Angstrom")
    parser.add_argument("--axis", choices=sorted(AXIS_TO_VECTOR), default="x", help="Shift axis")
    parser.add_argument("--pymol", default="pymol", help="PyMOL executable")
    parser.add_argument("--keep-temp", action="store_true", help="Keep generated .pdb and .pml files")
    return parser.parse_args()


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


def parse_target_residue(residue: str) -> Tuple[str, int]:
    match = TARGET_RE.match(residue.strip())
    if not match:
        raise ValueError(f"Residue must look like B0071 or C91, got: {residue!r}")
    chain, resid_s = match.groups()
    return chain.upper(), int(resid_s)


def parse_step2_pdb(path: Path) -> Step2Data:
    atoms_by_conformer: Dict[Tuple[str, str, int, str], List[AtomRecord]] = {}
    conf_codes: Dict[Tuple[str, str, int, str], str] = {}

    with path.open() as handle:
        for raw in handle:
            if not raw.startswith(("ATOM", "HETATM")):
                continue

            token = raw[21:30].strip()
            match = MCCE_RES_RE.match(token)
            if not match:
                continue

            chain, resid_s, state_id = match.groups()
            resid = int(resid_s)
            atom_name = raw[12:16].strip()
            resname = raw[17:20].strip()
            conf_code = parse_conf_code(raw)
            atom = AtomRecord(
                record=raw[0:6].strip() or "ATOM",
                atom_name=atom_name,
                resname=resname,
                chain=chain,
                resid=resid,
                state_id=state_id,
                conf_code=conf_code,
                x=float(raw[30:38]),
                y=float(raw[38:46]),
                z=float(raw[46:54]),
                element=infer_element(atom_name),
            )
            key = (resname, chain, resid, state_id)
            atoms_by_conformer.setdefault(key, []).append(atom)
            if conf_code:
                conf_codes.setdefault(key, conf_code)

    return Step2Data(atoms_by_conformer=atoms_by_conformer, conf_codes=conf_codes)


def parse_conf_code(raw: str) -> str:
    fields = raw[54:].split()
    if not fields:
        return ""
    tail = fields[-1]
    return tail[:2] if len(tail) >= 2 else ""


def parse_fort38_labels(path: Path) -> Dict[Tuple[str, int], List[str]]:
    labels: Dict[Tuple[str, int], List[str]] = {}
    if not path.exists():
        return labels

    with path.open() as handle:
        for raw in handle:
            fields = raw.split()
            if not fields:
                continue
            match = FORT38_LABEL_RE.match(fields[0])
            if not match:
                continue
            _resname, _code, chain, resid_s, _state_id = match.groups()
            labels.setdefault((chain, int(resid_s)), []).append(fields[0])

    return labels


def build_microstates(
    step2: Step2Data,
    chain: str,
    resid: int,
    fort38_labels: Optional[Dict[Tuple[str, int], List[str]]] = None,
) -> List[Microstate]:
    residue_conformers = {
        key: atoms
        for key, atoms in step2.atoms_by_conformer.items()
        if key[1] == chain and key[2] == resid
    }
    if not residue_conformers:
        raise ValueError(f"No conformers found for residue {chain}{resid:04d}")

    resnames = sorted({key[0] for key in residue_conformers})
    if len(resnames) != 1:
        joined = ", ".join(resnames)
        raise ValueError(f"Residue {chain}{resid:04d} has multiple residue names: {joined}")
    resname = resnames[0]

    backbone_key = (resname, chain, resid, "000")
    backbone = step2.atoms_by_conformer.get(backbone_key)
    if not backbone:
        raise ValueError(f"Missing backbone conformer _000 for {resname} {chain}{resid:04d}")

    sidechain_keys = sorted(
        [key for key in residue_conformers if key[3] != "000"],
        key=lambda key: int(key[3]),
    )
    if not sidechain_keys:
        raise ValueError(f"No sidechain microstates found for {resname} {chain}{resid:04d}")

    labels_by_state = labels_for_residue(step2, resname, chain, resid, sidechain_keys, fort38_labels)
    sidechain_keys = order_sidechain_keys_by_labels(sidechain_keys, chain, resid, fort38_labels)
    states: List[Microstate] = []
    for key in sidechain_keys:
        state_id = key[3]
        label = labels_by_state[state_id]
        atoms = tuple(backbone + step2.atoms_by_conformer[key])
        states.append(
            Microstate(
                label=label,
                object_name=safe_pymol_name(label),
                resname=resname,
                chain=chain,
                resid=resid,
                state_id=state_id,
                atoms=atoms,
            )
        )
    return states


def labels_for_residue(
    step2: Step2Data,
    resname: str,
    chain: str,
    resid: int,
    sidechain_keys: Sequence[Tuple[str, str, int, str]],
    fort38_labels: Optional[Dict[Tuple[str, int], List[str]]],
) -> Dict[str, str]:
    labels = (fort38_labels or {}).get((chain, resid), [])
    labels_by_state: Dict[str, str] = {}

    for label in labels:
        match = FORT38_LABEL_RE.match(label)
        if not match:
            continue
        label_resname, _code, _chain, _resid_s, state_id = match.groups()
        if label_resname == resname:
            labels_by_state[state_id] = label

    for key in sidechain_keys:
        state_id = key[3]
        if state_id in labels_by_state:
            continue
        code = step2.conf_codes.get(key, "")
        if not code:
            code = "??"
        labels_by_state[state_id] = f"{resname}{code}{chain}{resid:04d}_{state_id}"

    return labels_by_state


def order_sidechain_keys_by_labels(
    sidechain_keys: Sequence[Tuple[str, str, int, str]],
    chain: str,
    resid: int,
    fort38_labels: Optional[Dict[Tuple[str, int], List[str]]],
) -> List[Tuple[str, str, int, str]]:
    keys_by_state = {key[3]: key for key in sidechain_keys}
    ordered_keys: List[Tuple[str, str, int, str]] = []
    seen: set[str] = set()

    for label in (fort38_labels or {}).get((chain, resid), []):
        match = FORT38_LABEL_RE.match(label)
        if not match:
            continue
        state_id = match.group(5)
        key = keys_by_state.get(state_id)
        if key and state_id not in seen:
            ordered_keys.append(key)
            seen.add(state_id)

    for key in sidechain_keys:
        if key[3] not in seen:
            ordered_keys.append(key)

    return ordered_keys


def safe_pymol_name(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", label)


def default_out_prefix(first_microstate: Microstate) -> str:
    return f"{first_microstate.resname}_{first_microstate.chain}{first_microstate.resid:04d}_microstates"


def write_microstate_pdb(path: Path, microstate: Microstate) -> None:
    lines: List[str] = []
    for serial, atom in enumerate(microstate.atoms, start=1):
        lines.append(atom_to_pdb_line(atom, serial))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def atom_to_pdb_line(atom: AtomRecord, serial: int) -> str:
    atom_name = format_atom_name(atom.atom_name)
    element = atom.element.strip() or infer_element(atom.atom_name)
    return (
        f"{atom.record:<6}{serial:>5d} {atom_name} "
        f"{atom.resname:>3} {atom.chain:1}{atom.resid:>4d}    "
        f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {element:>2}"
    )


def render_pml(
    microstates: Sequence[Microstate],
    pdb_paths: Sequence[Path],
    pse_path: Path,
    spacing: float,
    axis: str,
) -> str:
    axis_vector = AXIS_TO_VECTOR[axis]
    lines = [
        "reinitialize",
        "set retain_order, 1",
        "set auto_zoom, off",
        "bg_color white",
    ]

    for index, (microstate, pdb_path) in enumerate(zip(microstates, pdb_paths)):
        offset = tuple(component * spacing * index for component in axis_vector)
        object_name = microstate.object_name
        label_name = f"{object_name}_label"
        label_pos = label_position(microstate, offset)
        lines.extend(
            [
                f"load {pdb_path.as_posix()}, {object_name}",
                f"translate [{offset[0]:.3f}, {offset[1]:.3f}, {offset[2]:.3f}], {object_name}",
                f"show sticks, {object_name}",
                f"hide lines, {object_name}",
                f"util.cbag {object_name}",
                (
                    f"pseudoatom {label_name}, pos=[{label_pos[0]:.3f}, {label_pos[1]:.3f}, "
                    f"{label_pos[2]:.3f}]"
                ),
                f'label {label_name}, "{microstate.label}"',
                f"set label_color, black, {label_name}",
            ]
        )

    all_objects = " ".join(ms.object_name for ms in microstates)
    lines.extend(
        [
            f"zoom {all_objects}" if all_objects else "zoom",
            f"save {pse_path.as_posix()}",
            "quit",
            "",
        ]
    )
    return "\n".join(lines)


def label_position(microstate: Microstate, offset: Tuple[float, float, float]) -> Tuple[float, float, float]:
    anchor = next((atom for atom in microstate.atoms if atom.atom_name == "CA"), microstate.atoms[0])
    return anchor.x + offset[0], anchor.y + offset[1], anchor.z + offset[2] + 2.0


def run_pymol(pymol: str, pml_path: Path) -> None:
    command = [pymol, "-cq", str(pml_path)]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"PyMOL executable not found: {pymol!r}") from exc
    except subprocess.CalledProcessError as exc:
        joined = " ".join(command)
        raise RuntimeError(f"PyMOL command failed with exit code {exc.returncode}: {joined}") from exc


def resolve_fort38(cli_value: Optional[str]) -> Optional[Path]:
    if cli_value:
        return Path(cli_value)
    default_path = Path("fort.38")
    return default_path if default_path.exists() else None


def export_residue_microstates(args: argparse.Namespace) -> Path:
    step2_path = Path(args.step2_pdb)
    chain, resid = parse_target_residue(args.residue)
    fort38_path = resolve_fort38(args.fort38)
    fort38_labels = parse_fort38_labels(fort38_path) if fort38_path else {}
    step2 = parse_step2_pdb(step2_path)
    microstates = build_microstates(step2, chain, resid, fort38_labels)

    out_prefix = args.out_prefix or default_out_prefix(microstates[0])
    pse_path = Path(f"{out_prefix}.pse")
    temp_root = Path.cwd() if args.keep_temp else Path(tempfile.mkdtemp(prefix=f"{out_prefix}_"))

    try:
        pdb_paths = []
        for microstate in microstates:
            pdb_path = temp_root / f"{out_prefix}_{microstate.object_name}.pdb"
            write_microstate_pdb(pdb_path, microstate)
            pdb_paths.append(pdb_path)

        pml_path = temp_root / f"{out_prefix}.pml"
        pml_path.write_text(render_pml(microstates, pdb_paths, pse_path, args.spacing, args.axis))
        run_pymol(args.pymol, pml_path)
    finally:
        if not args.keep_temp and temp_root.exists():
            shutil.rmtree(temp_root)

    return pse_path


def main() -> int:
    args = parse_args()
    try:
        pse_path = export_residue_microstates(args)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1
    print(f"Wrote {pse_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
