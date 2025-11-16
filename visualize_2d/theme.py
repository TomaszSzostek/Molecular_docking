"""
Color palettes and typography tokens for visualize_2d.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Palette:
    background: str = "#f5f7fb"
    panel: str = "#ffffff"
    accent: str = "#95b5f7"
    accent_secondary: str = "#f7b5c3"
    text_dark: str = "#23303d"
    text_muted: str = "#6b7684"
    outline: str = "#dce4f2"
    ligands: tuple[str, ...] = ("#ffb677", "#f089b9", "#7bdff2", "#90d26d")
    interactions: dict[str, str] = None  # type: ignore

    def __post_init__(self):
        if self.interactions is None:
            object.__setattr__(
                self,
                "interactions",
                {
                    "hb": "#ff9a8d",
                    "hp": "#ffd180",
                    "pc": "#9fa8ff",
                    "ps": "#8fd6ff",
                    "sb": "#a3ffd6",
                    "ha": "#ffc3a0",
                    "wb": "#b5e0ff",
                    "me": "#c7b0ff",
                },
            )


PALETTE = Palette()

TITLE_FONT = "Avenir Next, Inter, Helvetica, Arial, sans-serif"
BODY_FONT = "Inter, Helvetica, Arial, sans-serif"

