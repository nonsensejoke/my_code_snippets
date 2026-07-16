#!/usr/bin/env python3
import argparse
import math
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_PYMOL = "~/packages-installed/bin/pymol"
DEFAULT_OBABEL = "~/packages-installed/bin/obabel"
DEFAULT_XTB = "~/.miniconda3/envs/xtb-1/bin/xtb"


@dataclass
class Atom:
    x: float
    y: float
    z: float
    element: str
    charge_code: int
    line_index: int


@dataclass
class SdfFile:
    path: Path
    lines: List[str]
    atom_count: int
    bond_count: int
    atom_start: int
    atoms: List[Atom]

    def with_coordinates(self, coordinates: Sequence[Tuple[str, float, float, float]]) -> str:
        if len(coordinates) != self.atom_count:
            raise ValueError(
                f"coordinate count mismatch for {self.path}: "
                f"expected {self.atom_count}, got {len(coordinates)}"
            )

        updated = list(self.lines)
        for atom, coordinate in zip(self.atoms, coordinates):
            element, x, y, z = coordinate
            if element != atom.element:
                raise ValueError(
                    f"element mismatch at atom {atom.line_index + 1}: "
                    f"SDF has {atom.element}, coordinates have {element}"
                )
            updated[atom.line_index] = format_atom_line(updated[atom.line_index], x, y, z)
        return "\n".join(updated) + "\n"


@dataclass
class ToolPaths:
    pymol: str = DEFAULT_PYMOL
    obabel: str = DEFAULT_OBABEL
    xtb: str = DEFAULT_XTB


def expand_tool_paths(tools: ToolPaths) -> ToolPaths:
    return ToolPaths(
        pymol=str(Path(tools.pymol).expanduser()),
        obabel=str(Path(tools.obabel).expanduser()),
        xtb=str(Path(tools.xtb).expanduser()),
    )


@dataclass
class WorkflowResult:
    sdf0_output: Path
    sdf1_output: Path
    xyz_output: Optional[Path]
    log_output: Optional[Path]
    net_charge: int
    charge_sum: float
    workdir: Path


def parse_single_sdf(path: Path) -> SdfFile:
    text = path.read_text()
    if text.count("$$$$") != 1:
        raise ValueError(f"only single-molecule SDF files are supported: {path}")

    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError(f"not enough lines for SDF header: {path}")

    count_parts = lines[3].split()
    if len(count_parts) < 2 or "V2000" not in lines[3]:
        raise ValueError(f"only V2000 SDF counts lines are supported: {path}")

    atom_count = int(count_parts[0])
    bond_count = int(count_parts[1])
    atom_start = 4
    atom_end = atom_start + atom_count
    if len(lines) < atom_end:
        raise ValueError(f"SDF atom block is shorter than count line declares: {path}")

    atoms = [parse_atom_line(lines[index], index) for index in range(atom_start, atom_end)]
    return SdfFile(
        path=path,
        lines=lines,
        atom_count=atom_count,
        bond_count=bond_count,
        atom_start=atom_start,
        atoms=atoms,
    )


def parse_atom_line(line: str, line_index: int) -> Atom:
    parts = line.split()
    if len(parts) < 6:
        raise ValueError(f"cannot parse atom line {line_index + 1}: {line}")
    return Atom(
        x=float(parts[0]),
        y=float(parts[1]),
        z=float(parts[2]),
        element=parts[3],
        charge_code=int(parts[5]),
        line_index=line_index,
    )


def format_atom_line(original_line: str, x: float, y: float, z: float) -> str:
    suffix = original_line[30:] if len(original_line) >= 30 else " " + " ".join(original_line.split()[3:])
    return f"{x:10.4f}{y:10.4f}{z:10.4f}{suffix}"


def parse_partial_charge_report(report: str) -> Tuple[float, int]:
    charges: List[float] = []
    in_charge_block = False
    for raw_line in report.splitlines():
        line = raw_line.strip()
        if line == "ATOMIC CHARGES":
            in_charge_block = True
            continue
        if not in_charge_block:
            continue
        if not line:
            if charges:
                break
            continue
        parts = line.split()
        if len(parts) < 3:
            if charges:
                break
            continue
        try:
            charges.append(float(parts[-1]))
        except ValueError:
            if charges:
                break

    if not charges:
        raise ValueError("Open Babel report did not contain an ATOMIC CHARGES block")

    charge_sum = math.fsum(charges)
    return charge_sum, int(round(charge_sum))


def parse_xyz(path: Path) -> List[Tuple[str, float, float, float]]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"not enough lines for XYZ file: {path}")
    atom_count = int(lines[0].strip())
    atom_lines = lines[2 : 2 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(f"XYZ atom block is shorter than count line declares: {path}")

    coordinates: List[Tuple[str, float, float, float]] = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"cannot parse XYZ atom line: {line}")
        coordinates.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return coordinates


def write_xyz(path: Path, coordinates: Sequence[Tuple[str, float, float, float]], title: str) -> None:
    lines = [str(len(coordinates)), title]
    lines.extend(f"{element} {x:.8f} {y:.8f} {z:.8f}" for element, x, y, z in coordinates)
    path.write_text("\n".join(lines) + "\n")


def match_original_atoms_to_hydrogenated(
    original: SdfFile, hydrogenated: SdfFile, tolerance: float = 1.0e-3
) -> Dict[int, int]:
    available = set(range(hydrogenated.atom_count))
    mapping: Dict[int, int] = {}

    for original_index, atom in enumerate(original.atoms):
        matches = []
        for hydrogenated_index in available:
            candidate = hydrogenated.atoms[hydrogenated_index]
            if atom.element != candidate.element:
                continue
            distance = math.dist((atom.x, atom.y, atom.z), (candidate.x, candidate.y, candidate.z))
            if distance <= tolerance:
                matches.append(hydrogenated_index)
        if len(matches) != 1:
            raise ValueError(
                f"could not uniquely map original atom {original_index + 1} "
                f"({atom.element}) to hydrogenated SDF; matches={matches}"
            )
        mapping[original_index] = matches[0]
        available.remove(matches[0])

    return mapping


class CommandRunner:
    def add_hydrogens(self, input_sdf: Path, output_sdf: Path, tools: ToolPaths, workdir: Path) -> None:
        command = f"load {input_sdf}, mol; h_add mol; save {output_sdf}, mol; quit"
        self._run([tools.pymol, "-cq", "-d", command], cwd=workdir)
        if not output_sdf.exists():
            raise RuntimeError(f"PyMOL finished without creating {output_sdf}")

    def partial_charge_report(self, input_sdf: Path, tools: ToolPaths, workdir: Path) -> str:
        result = self._run(
            [tools.obabel, str(input_sdf), "-o", "report", "--partialcharge", "mmff94"],
            cwd=workdir,
        )
        return result.stdout

    def sdf_to_xyz(self, input_sdf: Path, output_xyz: Path, tools: ToolPaths, workdir: Path) -> None:
        self._run([tools.obabel, str(input_sdf), "-O", str(output_xyz)], cwd=workdir)

    def optimize_xtb(self, input_xyz: Path, net_charge: int, tools: ToolPaths, workdir: Path) -> Path:
        self._run(
            [tools.xtb, str(input_xyz), "--gfn", "2", "--opt", "--chrg", str(net_charge)],
            cwd=workdir,
        )
        output = workdir / "xtbopt.xyz"
        if not output.exists():
            raise RuntimeError(f"xTB finished without creating {output}")
        return output

    def _run(self, args: Sequence[str], cwd: Path) -> subprocess.CompletedProcess:
        return subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True)


class FakeCommandRunner(CommandRunner):
    def __init__(self, hydrogenated_sdf: Path, partial_charge_report: str, optimized_xyz: str):
        self.hydrogenated_sdf = hydrogenated_sdf
        self.partial_charge_report_text = partial_charge_report
        self.optimized_xyz_text = optimized_xyz
        self.hydrogenated_output: Path = Path()
        self.report_input: Path = Path()
        self.xtb_args: List[str] = []

    def add_hydrogens(self, input_sdf: Path, output_sdf: Path, tools: ToolPaths, workdir: Path) -> None:
        shutil.copyfile(self.hydrogenated_sdf, output_sdf)
        self.hydrogenated_output = output_sdf

    def partial_charge_report(self, input_sdf: Path, tools: ToolPaths, workdir: Path) -> str:
        self.report_input = input_sdf
        return self.partial_charge_report_text

    def sdf_to_xyz(self, input_sdf: Path, output_xyz: Path, tools: ToolPaths, workdir: Path) -> None:
        sdf = parse_single_sdf(input_sdf)
        write_xyz(output_xyz, [(atom.element, atom.x, atom.y, atom.z) for atom in sdf.atoms], input_sdf.stem)

    def optimize_xtb(self, input_xyz: Path, net_charge: int, tools: ToolPaths, workdir: Path) -> Path:
        self.xtb_args = [tools.xtb, str(input_xyz), "--gfn", "2", "--opt", "--chrg", str(net_charge)]
        output = workdir / "xtbopt.xyz"
        output.write_text(self.optimized_xyz_text)
        return output


def run_workflow(
    input_sdf: Path,
    out_dir: Optional[Path],
    tools: ToolPaths,
    runner: CommandRunner = None,
    keep_workdir: bool = False,
) -> WorkflowResult:
    runner = runner or CommandRunner()
    tools = expand_tool_paths(tools)
    input_sdf = input_sdf.resolve()
    out_dir = (out_dir or input_sdf.parent).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = input_sdf.stem
    workdir = Path(tempfile.mkdtemp(prefix=f"{stem}_", dir=str(out_dir)))

    try:
        original = parse_single_sdf(input_sdf)
        hydrogenated_path = workdir / f"{stem}_h.sdf"
        runner.add_hydrogens(input_sdf, hydrogenated_path, tools, workdir)

        hydrogenated = parse_single_sdf(hydrogenated_path)
        mapping = match_original_atoms_to_hydrogenated(original, hydrogenated)

        report = runner.partial_charge_report(hydrogenated_path, tools, workdir)
        charge_sum, net_charge = parse_partial_charge_report(report)

        xtb_input = workdir / f"{stem}_h.xyz"
        runner.sdf_to_xyz(hydrogenated_path, xtb_input, tools, workdir)
        optimized_xyz = runner.optimize_xtb(xtb_input, net_charge, tools, workdir)
        optimized_coordinates = parse_xyz(optimized_xyz)

        sdf1_text = hydrogenated.with_coordinates(optimized_coordinates)
        original_coordinates = [optimized_coordinates[mapping[index]] for index in range(original.atom_count)]
        sdf0_text = original.with_coordinates(original_coordinates)

        sdf0_output = out_dir / f"{stem}_xtbopt.sdf"
        sdf1_output = out_dir / f"{stem}_xtbopt-h.sdf"

        sdf0_output.write_text(sdf0_text)
        sdf1_output.write_text(sdf1_text)

        result = WorkflowResult(
            sdf0_output=sdf0_output,
            sdf1_output=sdf1_output,
            xyz_output=None,
            log_output=None,
            net_charge=net_charge,
            charge_sum=charge_sum,
            workdir=workdir,
        )
        return result
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add hydrogens, optimize with xTB GFN2, and update SDF coordinates.")
    parser.add_argument("input_sdf", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--pymol", default=DEFAULT_PYMOL)
    parser.add_argument("--obabel", default=DEFAULT_OBABEL)
    parser.add_argument("--xtb", default=DEFAULT_XTB)
    parser.add_argument("--keep-workdir", action="store_true")
    return parser


def main(argv: Sequence[str] = None) -> int:
    args = build_parser().parse_args(argv)
    tools = ToolPaths(pymol=args.pymol, obabel=args.obabel, xtb=args.xtb)
    try:
        result = run_workflow(args.input_sdf, args.out_dir, tools, keep_workdir=args.keep_workdir)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"net_charge={result.net_charge} partial_charge_sum={result.charge_sum:.6f}")
    print(f"sdf={result.sdf0_output}")
    print(f"sdf_h={result.sdf1_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
