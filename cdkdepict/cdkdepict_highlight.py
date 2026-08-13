#!/usr/bin/env python3

"""Download an atom-highlighted CDK depiction from the public CDK Depict service.

Highlighting works by tagging atoms with SMILES atom-map numbers and asking for
``annotate=colmap``.  The map number is not drawn -- it only selects a colour
from CDK's fixed palette (see PALETTE below).  Because the number travels with
the atom inside the SMILES, the highlight lands on the right atom regardless of
how CDK orders atoms internally.

    # built-in aspirin demo (stdlib only)
    python3 cdkdepict_highlight.py

    # highlight by 0-based atom index (needs RDKit)
    python3 cdkdepict_highlight.py \
      --smiles 'CC(=O)Oc1ccccc1C(=O)O' \
      --highlight '0,1' --highlight '3,4,5:5' \
      --output aspirin.png

    # already-mapped SMILES, no RDKit needed
    python3 cdkdepict_highlight.py --mapped-smiles '[CH3:1]C(=O)O'

    python3 cdkdepict_highlight.py --list-colors

Limitations of the public endpoint (both verified against the live service):

* Colours come from the fixed palette below; arbitrary #RRGGBB is not possible.
  A local JSON wrapper accepting arbitrary colours was tried and removed --
  stock cdkdepict exposes no such POST API, so that code path never worked.
* Map numbers above 20 are silently ignored: the service still returns a valid
  PNG, just with no highlight at all.  This script rejects them up front.
* ``annotate`` takes a single value, so one image cannot show both atom numbers
  and highlights.
"""

import argparse
import os
import sys
import urllib.error
import urllib.parse
import urllib.request


DEFAULT_BASE_URL = "https://www.simolecule.com/cdkdepict/depict"
DEFAULT_STYLE = "cow"
STYLES = ("cow", "cob", "cot", "bow", "bot", "wob", "nob")

# CDK's colmap palette, sampled from the live endpoint.  Map number 0 means
# "unmapped" (no highlight), so the palette effectively starts at 1.
PALETTE = {
    1: "#3CB44B",   2: "#FFE119",   3: "#0082C8",   4: "#F58231",
    5: "#911EB4",   6: "#46F0F0",   7: "#F032E6",   8: "#D2F53C",
    9: "#FABEBE",  10: "#008080",  11: "#E6BEFF",  12: "#AA6E28",
    13: "#FFFAC8", 14: "#800000",  15: "#AAFFC3",  16: "#808000",
    17: "#FFD8B1", 18: "#000080",  19: "#808080",  20: "#E3E3E3",
}
MAX_COLOR = max(PALETTE)

# Stdlib-only demo: map 4 paints atoms [0,1] orange, map 3 paints [3,4,5] blue.
DEFAULT_MAPPED_SMILES = "[CH3:4][C:4](=O)[O:3][c:3]1[cH:3][cH][cH][cH][c]1C(=O)O"


def parse_highlight(value):
    """Parse ``ATOM,ATOM`` or ``ATOM,ATOM:COLOR`` into ``(atoms, colour_or_None)``.

    COLOR is a palette number (1-%d), not a hex code -- the public endpoint
    cannot take arbitrary colours.  Groups without one get the next free number.
    """ % MAX_COLOR
    atoms_text, sep, color_text = value.partition(":")

    try:
        atoms = [int(item.strip()) for item in atoms_text.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "highlight must look like '0,1' or '0,1:3'"
        ) from exc
    if not atoms or any(atom < 0 for atom in atoms):
        raise argparse.ArgumentTypeError("atom indexes must be non-negative integers")

    if not sep:
        return atoms, None

    color_text = color_text.strip()
    if color_text.startswith("#"):
        raise argparse.ArgumentTypeError(
            "the public endpoint cannot take hex colours; use a palette number "
            "1-%d instead (see --list-colors)" % MAX_COLOR
        )
    try:
        color = int(color_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "colour must be a palette number 1-%d (see --list-colors)" % MAX_COLOR
        ) from exc
    if color not in PALETTE:
        raise argparse.ArgumentTypeError(
            "colour %d is out of range; the service silently drops highlights "
            "above %d, so pick 1-%d (see --list-colors)" % (color, MAX_COLOR, MAX_COLOR)
        )
    return atoms, color


def assign_colors(groups):
    """Fill in the colour of every group that did not ask for a specific one."""
    used = {color for _atoms, color in groups if color is not None}
    free = (n for n in sorted(PALETTE) if n not in used)
    out = []
    for atoms, color in groups:
        if color is None:
            color = next(free, None)
            if color is None:
                raise SystemExit("too many highlight groups (max %d)" % MAX_COLOR)
        out.append((atoms, color))
    return out


def build_mapped_smiles(smiles, groups):
    """Turn 0-based atom indexes into atom-map numbers.  Needs RDKit."""
    try:
        from rdkit import Chem, RDLogger
    except ImportError:
        raise SystemExit(
            "--smiles/--highlight needs RDKit (pip install rdkit); "
            "otherwise pass a pre-mapped SMILES via --mapped-smiles"
        )
    RDLogger.DisableLog("rdApp.*")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit("could not parse --smiles: %s" % smiles)
    count = mol.GetNumAtoms()
    for atoms, color in groups:
        for index in atoms:
            if index >= count:
                raise SystemExit(
                    "atom index %d is out of range (molecule has %d heavy atoms, "
                    "so valid indexes are 0-%d)" % (index, count, count - 1)
                )
            mol.GetAtomWithIdx(index).SetAtomMapNum(color)
    # canonical=False keeps the output order close to the input for readability;
    # correctness does not depend on it since the map rides along with the atom.
    return Chem.MolToSmiles(mol, canonical=False)


def validate_png(data):
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("response is not a PNG image")


def request_png(base_url, style, mapped_smiles, zoom, timeout):
    """GET a highlighted PNG using CDK's built-in colmap annotation."""
    query = urllib.parse.urlencode(
        {
            "smi": mapped_smiles,
            "w": "-1",
            "h": "-1",
            "abbr": "off",
            "hdisp": "S",
            "zoom": str(zoom),
            "annotate": "colmap",
            "r": "0",
        }
    )
    request_url = "%s/%s/png?%s" % (base_url.rstrip("/"), style, query)
    request = urllib.request.Request(
        request_url,
        headers={"Accept": "image/png", "User-Agent": "cdkdepict-highlight/3.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read()
        content_type = response.headers.get("Content-Type", "")
    validate_png(data)
    return data, content_type, request_url


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Download an atom-highlighted CDK PNG from the public service."
    )
    parser.add_argument("--smiles", help="SMILES to highlight by atom index (needs RDKit)")
    parser.add_argument(
        "--highlight",
        action="append",
        type=parse_highlight,
        metavar="ATOMS[:COLOR]",
        default=[],
        help="repeatable group of 0-based atom indexes, e.g. '0,1' or '3,4,5:5'; "
             "COLOR is a palette number 1-%d, auto-assigned when omitted" % MAX_COLOR,
    )
    parser.add_argument(
        "--mapped-smiles",
        help="SMILES that already carries atom-map numbers (skips RDKit)",
    )
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help="CDK Depict base URL")
    parser.add_argument("--style", default=DEFAULT_STYLE, choices=STYLES,
                        help="colour scheme (default: %(default)s)")
    parser.add_argument("--zoom", type=float, default=2.2, help="Depiction zoom")
    parser.add_argument("--timeout", type=float, default=60, help="Request timeout")
    parser.add_argument("-o", "--output", default="cdkdepict_highlight.png")
    parser.add_argument("--print-request", action="store_true")
    parser.add_argument("--list-colors", action="store_true",
                        help="print the palette (map number -> colour) and exit")
    args = parser.parse_args(argv)

    if args.list_colors:
        print("CDK colmap palette (atom map number -> colour):")
        for number in sorted(PALETTE):
            print("  %2d  %s" % (number, PALETTE[number]))
        print("map numbers above %d are silently ignored by the service" % MAX_COLOR)
        return 0

    if args.zoom <= 0 or args.timeout <= 0:
        parser.error("--zoom and --timeout must be greater than zero")
    if args.mapped_smiles and (args.smiles or args.highlight):
        parser.error("--mapped-smiles cannot be combined with --smiles/--highlight")

    if args.mapped_smiles:
        mapped_smiles = args.mapped_smiles
    elif args.smiles or args.highlight:
        if not args.smiles:
            parser.error("--highlight needs --smiles")
        groups = assign_colors(args.highlight) if args.highlight else []
        if not groups:
            parser.error("--smiles needs at least one --highlight group")
        mapped_smiles = build_mapped_smiles(args.smiles, groups)
    else:
        mapped_smiles = DEFAULT_MAPPED_SMILES      # stdlib-only demo

    try:
        data, content_type, request_url = request_png(
            args.base_url, args.style, mapped_smiles, args.zoom, args.timeout
        )
    except urllib.error.HTTPError as exc:
        detail = exc.read(1000).decode("utf-8", errors="replace")
        print("HTTP %d: %s" % (exc.code, detail), file=sys.stderr)
        return 1
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 1

    if args.print_request:
        print(request_url)

    output = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    with open(output, "wb") as file:
        file.write(data)

    print("smiles: %s" % mapped_smiles)
    print("saved: %s" % output)
    print("bytes: %d" % len(data))
    print("content-type: %s" % content_type)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
