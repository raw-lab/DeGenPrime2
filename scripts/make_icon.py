#!/usr/bin/env python3
"""Regenerate DeGenPrime desktop icons (PNG + ICO + ICNS).

Run: `python3 scripts/make_icon.py` from the project root.
Requires Pillow (pip install pillow).
"""
from __future__ import annotations

import os
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    raise SystemExit("Pillow not installed. pip install pillow")

OUT_DIR = Path(__file__).resolve().parent.parent / "assets"


def make_icon(size: int) -> Image.Image:
    """Draw a DeGenPrime icon at the given pixel size.

    Style: layered blue circles, "DP" monogram in white, DNA-helix-style
    dots along the bottom.
    """
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)

    margin = max(1, size // 14)
    rings = [(43, 124, 224, 255), (60, 145, 230, 255), (90, 170, 240, 255)]
    for i, color in enumerate(rings):
        off = margin + i * max(1, size // 64)
        if size - off > off:
            d.ellipse([off, off, size - off, size - off], fill=color)

    if size >= 32:
        try:
            font = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                size // 2,
            )
        except Exception:
            font = ImageFont.load_default()
        text = "DP"
        bbox = d.textbbox((0, 0), text, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        d.text(
            ((size - tw) / 2 - bbox[0], (size - th) / 2 - bbox[1] - size // 20),
            text, fill="white", font=font,
        )

    if size >= 64:
        cx, cy = size // 2, int(size * 0.78)
        rad = max(1, size // 30)
        for i in range(-3, 4):
            x = cx + i * size // 16
            y = cy + int((size // 50) * (-1 if i % 2 else 1))
            d.ellipse([x - rad, y - rad, x + rad, y + rad], fill="white")
    return img


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Single-size PNG for Linux .desktop
    make_icon(256).save(OUT_DIR / "icon.png", "PNG")
    print(f"Wrote {OUT_DIR / 'icon.png'}")

    # Multi-size ICO for Windows
    sizes = [32, 48, 64, 128, 256]
    imgs = [make_icon(s) for s in sizes]
    imgs[-1].save(
        OUT_DIR / "icon.ico", format="ICO",
        sizes=[(s, s) for s in sizes],
    )
    print(f"Wrote {OUT_DIR / 'icon.ico'}")

    # Single-size ICNS for macOS
    try:
        make_icon(512).save(OUT_DIR / "icon.icns", format="ICNS")
        print(f"Wrote {OUT_DIR / 'icon.icns'}")
    except Exception as e:
        print(f"Skip .icns ({e}); macOS app will use a generic icon.")


if __name__ == "__main__":
    main()
