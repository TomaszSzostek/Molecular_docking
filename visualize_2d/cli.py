"""
CLI entry points for the pastel 2D interaction diagrams.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .data import ComplexAssets, build_complex_assets
from .drawer import render_board


def _cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="python -m visualize_2d",
        description="Render pastel 2D interaction boards for docking complexes.",
    )
    parser.add_argument(
        "--complex",
        type=Path,
        required=True,
        help="Path to the complex folder produced by run_pipeline/analyze files_for_visualization.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output image (.png or .svg). Default: <complex>/interaction_board.png",
    )
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        metavar=("W", "H"),
        default=(1600, 1300),
        help="Canvas size in pixels.",
    )
    parser.add_argument(
        "--format",
        choices=("png", "svg"),
        default=None,
        help="Override format (otherwise derived from --out suffix).",
    )
    return parser.parse_args()


def render_single(
    complex_dir: Path,
    out_path: Path | None = None,
    size: tuple[int, int] = (2200, 1800),
    fmt: str | None = None,
) -> Path:
    assets: ComplexAssets = build_complex_assets(complex_dir)
    out = out_path or (complex_dir / "interaction_board.png")
    fmt = fmt or out.suffix.lstrip(".") or "png"
    render_board(assets, out, size=size, fmt=fmt)
    return out


def main() -> None:
    args = _cli()
    render_single(args.complex, args.out, tuple(args.size), args.format)


if __name__ == "__main__":
    main()

