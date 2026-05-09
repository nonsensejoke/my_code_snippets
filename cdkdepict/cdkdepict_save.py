#!/usr/bin/env python3

"""
 用法：

  python3 scripts/cdkdepict_save.py "c1ccccc1(C) sample_molecule"

  指定格式和放大倍数：

  python3 scripts/cdkdepict_save.py "c1ccccc1(C) sample_molecule" --format png
  --zoom 1.3 -o sample.png
  python3 scripts/cdkdepict_save.py "c1ccccc1(C) sample_molecule" --format svg
  --zoom 2.0 -o sample.svg
  python3 scripts/cdkdepict_save.py "c1ccccc1(C) sample_molecule" --format pdf
  --zoom 1.0 -o sample.pdf
"""

import argparse
import os
import re
import sys
import urllib.parse
import urllib.request


BASE_URL = "https://www.simolecule.com/cdkdepict/depict"
FORMATS = {"png", "svg", "pdf"}
STYLES = {"cow", "cob", "cot", "bow", "bot", "wob", "nob"}


def sanitize_filename(value):
    cleaned = re.sub(r'[\\/:*?"<>|]+', "_", value.strip())
    cleaned = re.sub(r"\s+", "_", cleaned)
    cleaned = cleaned.strip("._")
    return cleaned or "molecule"


def title_from_smiles_line(smiles_line):
    parts = smiles_line.strip().split(maxsplit=1)
    if len(parts) == 2 and parts[1].strip():
        return parts[1].strip()
    return "molecule"


def build_url(smiles_line, image_format, zoom, style, abbr):
    query = urllib.parse.urlencode(
        {
            "smi": smiles_line,
            "w": "-1",
            "h": "-1",
            "abbr": abbr,
            "hdisp": "S",
            "zoom": str(zoom),
            "annotate": "none",
            "r": "0",
        }
    )
    return f"{BASE_URL}/{style}/{image_format}?{query}"


def resolve_output_path(output, title, image_format):
    filename = f"{sanitize_filename(title)}.{image_format}"
    if not output:
        return os.path.abspath(filename)

    output = os.path.abspath(output)
    if os.path.isdir(output) or output.endswith(os.sep):
        return os.path.join(output, filename)

    ext = os.path.splitext(output)[1].lower()
    if ext != f".{image_format}":
        os.makedirs(output, exist_ok=True)
        return os.path.join(output, filename)

    return output


def validate_payload(payload, image_format):
    if image_format == "png" and not payload.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("response is not a PNG file")
    if image_format == "pdf" and not payload.startswith(b"%PDF"):
        raise ValueError("response is not a PDF file")
    if image_format == "svg":
        head = payload[:300].lstrip()
        if not (head.startswith(b"<svg") or head.startswith(b"<?xml")):
            raise ValueError("response is not an SVG file")


def download(url, output_path, image_format):
    req = urllib.request.Request(url, headers={"User-Agent": "cdkdepict-save/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = resp.read()
        content_type = resp.headers.get("Content-Type", "")

    validate_payload(payload, image_format)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "wb") as f:
        f.write(payload)

    return len(payload), content_type


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Save a CDK Depict image without opening a browser.",
    )
    parser.add_argument(
        "smiles",
        help='SMILES line, optionally followed by molecule name, e.g. "c1ccccc1(C) sample_molecule"',
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=sorted(FORMATS),
        default="png",
        help="Output image format. Default: png",
    )
    parser.add_argument(
        "-z",
        "--zoom",
        type=float,
        default=1.3,
        help="Depiction zoom factor. Default: 1.3",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file or directory. Default: ./<molecule-title>.<format>",
    )
    parser.add_argument(
        "--abbr",
        choices=["on", "reagents", "groups", "off"],
        default="off",
        help='Abbreviation mode. Default: off ("Do Not Abbreviate")',
    )
    parser.add_argument(
        "--style",
        choices=sorted(STYLES),
        default="bot",
        help='CDK Depict style path. Default: bot ("Black on Clear")',
    )
    parser.add_argument(
        "--print-url",
        action="store_true",
        help="Print the request URL before downloading.",
    )

    args = parser.parse_args(argv)
    smiles_line = args.smiles.strip()
    if not smiles_line:
        parser.error("smiles cannot be empty")
    if args.zoom <= 0:
        parser.error("--zoom must be greater than 0")

    title = title_from_smiles_line(smiles_line)
    url = build_url(smiles_line, args.format, args.zoom, args.style, args.abbr)
    output_path = resolve_output_path(args.output, title, args.format)

    if args.print_url:
        print(url)

    try:
        size, content_type = download(url, output_path, args.format)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"saved: {output_path}")
    print(f"bytes: {size}")
    print(f"content-type: {content_type}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
