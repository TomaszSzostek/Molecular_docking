"""
Rendering logic for pastel 2D interaction boards.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
from typing import Tuple, Any
import math

from PIL import Image, ImageDraw, ImageFont, ImageFilter
from rdkit.Chem.Draw import rdMolDraw2D

from .data import ComplexAssets, INTERACTION_KIND_MAP, Interaction
from .theme import PALETTE


def _load_font(size: int, weight: str = "regular") -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return a system font with best-effort fallback."""
    families = {
        "regular": [
            "/System/Library/Fonts/SFNSRounded.ttf",
            "/System/Library/Fonts/Supplemental/Avenir.ttc",
            "/Library/Fonts/Arial.ttf",
        ],
        "bold": [
            "/Library/Fonts/Arial Bold.ttf",
            "/System/Library/Fonts/Supplemental/Avenir Next.ttc",
        ],
    }
    preferred = families.get(weight, []) + families.get("regular", [])
    for path in preferred:
        try:
            return ImageFont.truetype(path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def _render_ligand(mol, size: Tuple[int, int], highlight_atoms: set[int], highlight_colors: dict[int, tuple[float, float, float]]) -> tuple[Image.Image, dict[int, Tuple[float, float]]]:
    """Render the ligand to a Cairo surface and return image plus 2D atom coordinates."""
    w, h = size
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    options = drawer.drawOptions()
    options.fixedBondLength = 30
    options.padding = 0.05
    options.bondLineWidth = 5.0
    setattr(options, "atomLabelFontSize", 2.0)
    setattr(options, "atomLabelFontFace", "Arial")
    options.useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=list(highlight_atoms),
        highlightAtomColors=highlight_colors,
        highlightAtomRadii={idx: 0.6 for idx in highlight_atoms},

    )
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()
    img = Image.open(BytesIO(png_bytes)).convert("RGBA")
    coords = {idx: drawer.GetDrawCoords(idx) for idx in range(mol.GetNumAtoms())}
    coord_map = {idx: (pt.x, pt.y) for idx, pt in coords.items()}
    return img, coord_map


def _rounded_panel(draw: ImageDraw.ImageDraw, xy, radius: int, fill: str, outline: str | None = None):
    """Draw a rounded rectangle with optional outline."""
    x0, y0, x1, y1 = xy
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=2 if outline else 0)


def _hex_to_rgba(hex_color: str, alpha: int = 255) -> tuple[int, int, int, int]:
    """Convert #RRGGBB to an (r, g, b, a) tuple."""
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return (r, g, b, alpha)


def _draw_dashed_line(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 4, dash: int = 10):
    """Draw a dashed line between two points."""
    x0, y0 = start
    x1, y1 = end
    total = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
    if total == 0:
        return
    steps = int(total // dash)
    for i in range(0, steps, 2):
        t0 = i * dash / total
        t1 = min((i + 1) * dash / total, 1)
        sx = int(x0 + (x1 - x0) * t0)
        sy = int(y0 + (y1 - y0) * t0)
        ex = int(x0 + (x1 - x0) * t1)
        ey = int(y0 + (y1 - y0) * t1)
        draw.line([(sx, sy), (ex, ey)], fill=color, width=width)


def render_board(assets: ComplexAssets, out_path: Path, size: tuple[int, int] = (2200, 1800), fmt: str = "png") -> None:
    """Render a publication-style 2D interaction board for a docking complex."""
    width, height = size
    board = Image.new("RGBA", (width, height), PALETTE.background)
    draw = ImageDraw.Draw(board)

    # Panel with soft shadow
    panel = (80, 60, width - 80, height - 80)
    shadow = Image.new("RGBA", board.size, (0, 0, 0, 0))
    shadow_draw = ImageDraw.Draw(shadow)
    _rounded_panel(shadow_draw, (panel[0] + 8, panel[1] + 12, panel[2] + 8, panel[3] + 12), 40, "#00000020")
    shadow = shadow.filter(ImageFilter.GaussianBlur(18))
    board = Image.alpha_composite(board, shadow)
    draw = ImageDraw.Draw(board)
    _rounded_panel(draw, panel, 40, PALETTE.panel, PALETTE.outline)

    base_canvas = (2200, 1800)
    layout_scale = min(width / base_canvas[0], height / base_canvas[1])
    layout_scale = max(layout_scale, 0.5)

    def scale_metric(value: float, min_px: float = 1.0) -> int:
        return int(round(max(value * layout_scale, min_px)))

    # Ligand drawing
    ligand_box_width = scale_metric(1080, 640)
    ligand_box_height = scale_metric(780, 520)
    ligand_box_left = (width - ligand_box_width) // 2
    ligand_box_top = panel[1] + scale_metric(100, 60)
    ligand_box = (ligand_box_left, ligand_box_top, ligand_box_left + ligand_box_width, ligand_box_top + ligand_box_height)
    _rounded_panel(draw, ligand_box, 36, "#f9fbff")
    involved_atoms: set[int] = set()
    atom_colors: dict[int, tuple[float, float, float]] = {}
    for interaction in assets.interactions:
        for lig_atom in interaction.ligand_atoms:
            idx = assets.ligand_atom_map.get(lig_atom.strip())
            if idx is None:
                continue
            involved_atoms.add(idx)
            hex_color = PALETTE.interactions.get(interaction.kind, PALETTE.accent)
            rgb = tuple(int(hex_color[i : i + 2], 16) / 255 for i in (1, 3, 5))
            atom_colors[idx] = rgb

    ligand_img, atom_coords = _render_ligand(
        assets.ligand,
        (ligand_box[2] - ligand_box[0], ligand_box[3] - ligand_box[1]),
        involved_atoms,
        atom_colors,
    )
    board.paste(ligand_img, ligand_box[:2], ligand_img)

    # Titles (draw AFTER overlays so they remain visible)
    title_font = _load_font(56)
    subtitle_font = _load_font(34)
    title = f"{assets.ligand_id} ↦ {assets.receptor_id}"
    title_w = draw.textlength(title, font=title_font)
    draw.text(((width - title_w) // 2, panel[1] + 40), title, fill=PALETTE.text_dark, font=title_font)
    subtitle = "Detailed 2D interaction map"
    subtitle_w = draw.textlength(subtitle, font=subtitle_font)
    draw.text(((width - subtitle_w) // 2, panel[1] + 120), subtitle, fill="#4b4b4b", font=subtitle_font)

    # Interaction ribbons
    chip_font = _load_font(28)
    badge_draw = ImageDraw.Draw(board, "RGBA")
    center = ((ligand_box[0] + ligand_box[2]) // 2, (ligand_box[1] + ligand_box[3]) // 2)
    base_badge_height = scale_metric(68, 44)
    legend_forbidden_top = panel[3] - scale_metric(230, 150)
    badge_padding = scale_metric(40, 18)
    panel_guard = scale_metric(40, 24)
    vertical_guard = scale_metric(110, 70)
    ligand_guard = scale_metric(70, 42)
    radial_margin = scale_metric(50, 24)
    placed_boxes: list[tuple[float, float, float, float]] = []
    sector_count = 3
    sector_width = 2 * math.pi / sector_count
    sector_margin = sector_width * 0.08

    def _normalize_angle(angle: float) -> float:
        norm = angle % (2 * math.pi)
        return norm if norm >= 0 else norm + 2 * math.pi

    def _sector_center(idx: int) -> float:
        return (idx + 0.5) * sector_width

    def _angle_distance(a: float, b: float) -> float:
        diff = abs(a - b) % (2 * math.pi)
        return min(diff, 2 * math.pi - diff)

    def _sector_index(angle_norm: float) -> int:
        return int(angle_norm // sector_width) % sector_count

    def _angle_in_sector(theta_norm: float, sector_idx: int, loosen: int = 0) -> bool:
        allowance = (sector_width / 2) + sector_margin + loosen * sector_margin
        return _angle_distance(theta_norm, _sector_center(sector_idx)) <= allowance

    def _expanded_box(box):
        x0, y0, x1, y1 = box
        return (
            x0 - badge_padding,
            y0 - badge_padding,
            x1 + badge_padding,
            y1 + badge_padding,
        )

    def _boxes_overlap(box) -> bool:
        ax0, ay0, ax1, ay1 = _expanded_box(box)
        for bx0, by0, bx1, by1 in placed_boxes:
            ex0, ey0, ex1, ey1 = _expanded_box((bx0, by0, bx1, by1))
            if not (ax1 < ex0 or ax0 > ex1 or ay1 < ey0 or ay0 > ey1):
                return True
        return False

    def _within_bounds(box):
        x0, y0, x1, y1 = box
        return (
            x0 >= panel[0] + panel_guard
            and x1 <= panel[2] - panel_guard
            and y0 >= panel[1] + vertical_guard
            and y1 <= legend_forbidden_top
            and not (
                x1 >= ligand_box[0] - ligand_guard
                and x0 <= ligand_box[2] + ligand_guard
                and y1 >= ligand_box[1] - ligand_guard
                and y0 <= ligand_box[3] + ligand_guard
            )
        )

    def _clamp_box(box):
        x0, y0, x1, y1 = box
        w = x1 - x0
        h = y1 - y0
        x0 = max(panel[0] + panel_guard, min(panel[2] - panel_guard - w, x0))
        x1 = x0 + w
        y0 = max(panel[1] + vertical_guard, min(legend_forbidden_top - h, y0))
        y1 = y0 + h
        return [x0, y0, x1, y1]

    def _place_box(angle: float, angle_norm: float, badge_width: float, preferred_radius: float, sector_idx: int) -> tuple[float, float, float, float]:
        base_radius = scale_metric(55, 28)
        ring_step = scale_metric(35, 20)
        slots_per_ring = max(10, int(round(10 * layout_scale)))
        max_rings = 8
        preferred_ring = int(max(0, min(max_rings - 1, round((preferred_radius - base_radius) / ring_step))))
        ring_order = list(range(max_rings))
        ring_order.sort(key=lambda r: abs(r - preferred_ring))
        for ring in ring_order:
            radius = base_radius + ring * ring_step
            slots = slots_per_ring + ring * 4
            slot_angles = [(2 * math.pi * s) / slots for s in range(slots)]
            slot_angles.sort(key=lambda a: _angle_distance(_normalize_angle(a), angle_norm))
            for loosen in range(sector_count):
                placed = False
                for theta in slot_angles:
                    theta_norm = _normalize_angle(theta)
                    if not _angle_in_sector(theta_norm, sector_idx, loosen=loosen):
                        continue
                    cx = center[0] + radius * math.cos(theta)
                    cy = center[1] + radius * math.sin(theta)
                    box = _clamp_box(
                        (
                            cx - badge_width / 2,
                            cy - base_badge_height / 2,
                            cx + badge_width / 2,
                            cy + base_badge_height / 2,
                        )
                    )
                    if _within_bounds(box) and not _boxes_overlap(box):
                        return box
                    placed = True
                if placed:
                    break
        # fallback spiral if all rings busy
        radius = base_radius + max_rings * ring_step
        for attempt in range(100):
            theta = _sector_center(sector_idx) + ((-1) ** attempt) * (attempt / 12.0)
            cx = center[0] + radius * math.cos(theta)
            cy = center[1] + radius * math.sin(theta)
            box = _clamp_box(
                (
                    cx - badge_width / 2,
                    cy - base_badge_height / 2,
                    cx + badge_width / 2,
                    cy + base_badge_height / 2,
                )
            )
            if _within_bounds(box) and not _boxes_overlap(box):
                return box
            radius += scale_metric(20, 10)
        return _clamp_box(
            (
                center[0] + radius,
                center[1] + radius,
                center[0] + radius + badge_width,
                center[1] + radius + base_badge_height,
            )
        )

    def _mean_coords(atom_indices: list[int]):
        xs, ys = [], []
        for idx in atom_indices:
            if idx in atom_coords:
                xs.append(atom_coords[idx][0] + ligand_box[0])
                ys.append(atom_coords[idx][1] + ligand_box[1])
        if not xs:
            return None
        return (int(sum(xs) / len(xs)), int(sum(ys) / len(ys)))

    # group pi interactions by protein atom to avoid duplicates
    interactions: list[Interaction] = []
    seen_pi = {}
    for inter in assets.interactions:
        if inter.kind in {"pc", "ps"}:
            key = (inter.kind, inter.protein_atom, inter.protein_label)
            entry = seen_pi.setdefault(key, inter)
            if entry is not inter:
                entry.ligand_atoms = list(set(entry.ligand_atoms + inter.ligand_atoms))
            else:
                interactions.append(inter)
        else:
            interactions.append(inter)

    for interaction in interactions:
        atom_indices = []
        atom_points: list[tuple[int, float, float]] = []
        for lig_atom in interaction.ligand_atoms:
            idx = assets.ligand_atom_map.get(str(lig_atom).strip())
            if idx is None:
                continue
            atom_indices.append(idx)
            if idx in atom_coords:
                abs_x = atom_coords[idx][0] + ligand_box[0]
                abs_y = atom_coords[idx][1] + ligand_box[1]
                atom_points.append((idx, abs_x, abs_y))
        if not atom_points:
            continue
        centroid = (
            sum(p[1] for p in atom_points) / len(atom_points),
            sum(p[2] for p in atom_points) / len(atom_points),
        )
        anchor_idx, anchor_x, anchor_y = max(
            atom_points,
            key=lambda p: (p[1] - center[0]) ** 2 + (p[2] - center[1]) ** 2,
        )
        dx = anchor_x - center[0]
        dy = anchor_y - center[1]
        if dx == 0 and dy == 0:
            dx = 1
        norm = (dx ** 2 + dy ** 2) ** 0.5 or 1
        direction = (dx / norm, dy / norm)
        color = PALETTE.interactions.get(interaction.kind, PALETTE.accent)
        label_font = _load_font(24, weight="bold")
        residue = interaction.protein_label or interaction.protein_atom
        txt_w = badge_draw.textlength(residue, font=label_font)
        badge_width = max(int(txt_w) + 36, 72)
        angle = math.atan2(direction[1], direction[0])
        angle_norm = _normalize_angle(angle)
        sector_idx = _sector_index(angle_norm)
        preferred_radius = max(
            math.hypot(dx, dy) + ligand_guard + radial_margin,
            scale_metric(60, 30),
        )
        bx0, by0, bx1, by1 = _place_box(angle, angle_norm, badge_width, preferred_radius, sector_idx)
        placed_boxes.append((bx0, by0, bx1, by1))
        badge_center = ((bx0 + bx1) / 2, (by0 + by1) / 2)
        # Draw dashed lines to all interaction atoms (for pi-interactions with multiple atoms)
        for idx, atom_x, atom_y in atom_points:
            _draw_dashed_line(
                badge_draw,
                (int(atom_x), int(atom_y)),
                badge_center,
                color,
                width=3,
                dash=10,
            )
        badge_draw.rounded_rectangle((bx0, by0, bx1, by1), radius=14, fill=_hex_to_rgba(color, 205))
        badge_draw.rounded_rectangle((bx0 + 4, by0 + 4, bx1 - 4, by1 - 4), radius=10, outline=_hex_to_rgba("#FFFFFF", 150), width=2)
        text_origin = (badge_center[0] - txt_w / 2, badge_center[1] - label_font.size / 2 + 2)
        badge_draw.text(text_origin, residue, fill=PALETTE.text_dark, font=label_font)
        for idx, atom_x, atom_y in atom_points:
            px = int(atom_x)
            py = int(atom_y)
            badge_draw.ellipse((px - 8, py - 8, px + 8, py + 8), fill=_hex_to_rgba(color, 255))

    # Legend footer two rows
    legend_font = _load_font(28, weight="bold")
    legend_y = panel[3] - 150
    legend_title = "Legend"
    draw.text(((width - draw.textlength(legend_title, font=legend_font)) // 2, legend_y - 10), legend_title, fill=PALETTE.text_dark, font=legend_font)
    items = list(PALETTE.interactions.items())
    rows = 2
    cols = math.ceil(len(items) / rows)
    col_width = 300
    total_width = cols * col_width
    x_start = (width - total_width) // 2
    for idx, (kind, color) in enumerate(items):
        row = idx // cols
        col = idx % cols
        x = x_start + col * col_width
        y_offset = legend_y + 38 + row * 58
        human = INTERACTION_KIND_MAP.get(kind, (kind, kind))[0]
        box_size = 36
        draw.rounded_rectangle((x, y_offset, x + box_size, y_offset + box_size), radius=8, fill=color)
        draw.text((x + box_size + 10, y_offset + 4), human, fill=PALETTE.text_dark, font=_load_font(22))

    board = board.convert("RGB")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dpi = (800, 800)
    if fmt.lower() == "png":
        board.save(out_path, format="PNG", dpi=dpi)
    else:
        board.save(out_path, format=fmt.upper())

