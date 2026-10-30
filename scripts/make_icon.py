#!/usr/bin/env python3
"""Build every DeGenPrime icon and logo file from the master logo.

    python3 scripts/make_icon.py          (needs Pillow: pip install pillow)

Input:  assets/logo-source.png   the DeGenPrime robot logo with its wordmark
Output:
  assets/logo.png                       the logo, trimmed, 600 px wide (README, documents)
  assets/icon.png                       512 px app icon
  assets/icon.ico                       Windows icon, 16-256 px
  assets/icon.icns                      macOS icon, 16-1024 px
  assets/icons/degenprime-<n>.png       Linux icon theme sizes, 16-512 px
  web/assets/logo.webp                  the logo for the app's About page
  web/assets/mark.png                   the badge in the app's header
  web/assets/favicon.ico                browser tab icon
  web/assets/icon-192.png, icon-512.png install icons (web app manifest)
  web/assets/apple-touch-icon.png       iOS home-screen icon (opaque)

The app icon is the robot on a dark rounded square. Below 64 px only its head
is shown, which stays recognisable at 16 px where the whole bust does not.
After running this, scripts/make_windows_resource.sh embeds icon.ico in the
Windows executable.
"""
from __future__ import annotations

import io
import struct
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFilter
except ImportError:  # pragma: no cover
    raise SystemExit("Pillow is not installed: pip install pillow")

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "assets" / "logo-source.png"

# Regions of the 1254 x 1254 master: the robot above the wordmark, and its head.
BUST = (61, 32, 1217, 968)
HEAD = (440, 20, 1010, 590)
SMALL = 48  # sizes up to this use the head


def load() -> Image.Image:
    im = Image.open(SRC).convert("RGBA")
    if im.size != (1254, 1254):
        im = im.resize((1254, 1254), Image.LANCZOS)
    r, g, b, a = im.split()
    # The master's "opaque" pixels are 251-254; make them fully opaque.
    a = a.point(lambda v: 255 if v >= 245 else v)
    return Image.merge("RGBA", (r, g, b, a))


def gradient(size: int) -> Image.Image:
    top, bottom = (29, 47, 74), (10, 18, 32)
    col = Image.new("RGBA", (1, size))
    for y in range(size):
        t = y / max(1, size - 1)
        col.putpixel((0, y), tuple(round(top[i] * (1 - t) + bottom[i] * t) for i in range(3)) + (255,))
    bg = col.resize((size, size))
    # A soft teal glow behind the head, echoing the visor.
    glow = Image.new("L", (size, size), 0)
    ImageDraw.Draw(glow).ellipse([size * 0.18, size * 0.02, size * 0.82, size * 0.62], fill=110)
    glow = glow.filter(ImageFilter.GaussianBlur(size * 0.12))
    teal = Image.new("RGBA", (size, size), (0, 190, 220, 255))
    return Image.composite(Image.blend(bg, teal, 0.35), bg, glow)


def squircle(size: int) -> Image.Image:
    big = size * 4
    m = Image.new("L", (big, big), 0)
    ImageDraw.Draw(m).rounded_rectangle([0, 0, big - 1, big - 1], radius=int(big * 0.225), fill=255)
    return m.resize((size, size), Image.LANCZOS)


def badge(size: int, art: Image.Image, fill: float, rounded: bool = True) -> Image.Image:
    """The art on the dark background, bottom-aligned, as a rounded square."""
    canvas = gradient(size)
    w = max(1, round(size * fill))
    h = max(1, round(art.height * w / art.width))
    a = art.resize((w, h), Image.LANCZOS)
    canvas.alpha_composite(a, ((size - w) // 2, size - h))
    if not rounded:
        return canvas
    if size >= 64:
        # A faint light edge so the icon reads on dark docks and taskbars.
        edge = Image.new("L", (size, size), 0)
        r = int(size * 0.225)
        ImageDraw.Draw(edge).rounded_rectangle([0, 0, size - 1, size - 1], radius=r, outline=40, width=max(1, size // 128))
        canvas = Image.composite(Image.new("RGBA", (size, size), (255, 255, 255, 255)), canvas, edge)
    out = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    out.paste(canvas, (0, 0), squircle(size))
    return out


def icon(size: int, master: Image.Image) -> Image.Image:
    if size <= SMALL:
        return badge(size, master.crop(HEAD), 0.98)
    return badge(size, master.crop(BUST), 1.0)


def trimmed_logo(master: Image.Image, width: int | None = None) -> Image.Image:
    box = master.getchannel("A").getbbox()
    pad = 12
    box = (max(0, box[0] - pad), max(0, box[1] - pad), min(master.width, box[2] + pad), min(master.height, box[3] + pad))
    logo = master.crop(box)
    if width:
        logo = logo.resize((width, round(logo.height * width / logo.width)), Image.LANCZOS)
    return logo


def bmp_frame(im: Image.Image) -> bytes:
    """A 32-bit BGRA icon bitmap with its 1-bit transparency mask."""
    w, h = im.size
    rgba = im.convert("RGBA")
    px, alpha = rgba.tobytes(), rgba.getchannel("A").tobytes()
    xor, mask = bytearray(), bytearray()
    stride = (w + 31) // 32 * 4
    for y in range(h - 1, -1, -1):  # rows are stored bottom-up
        row = bytearray(px[y * w * 4 : (y + 1) * w * 4])
        row[0::4], row[2::4] = row[2::4], row[0::4]  # RGBA -> BGRA
        xor += row
        bits = bytearray(stride)
        for x in range(w):
            if alpha[y * w + x] < 128:
                bits[x // 8] |= 0x80 >> (x % 8)
        mask += bits
    header = struct.pack("<IiiHHIIiiII", 40, w, 2 * h, 1, 32, 0, len(xor) + len(mask), 0, 0, 0, 0)
    return header + bytes(xor) + bytes(mask)


def save_ico(images: list[Image.Image], path: Path) -> None:
    """A Windows icon: bitmap frames below 256 px, which every Windows API and
    resource compiler reads, and a PNG frame at 256 px, as Windows expects."""
    frames = []
    for im in images:
        if im.width >= 256:
            buf = io.BytesIO()
            im.save(buf, "PNG", optimize=True)
            frames.append((im.width, im.height, buf.getvalue()))
        else:
            frames.append((im.width, im.height, bmp_frame(im)))
    out = bytearray(struct.pack("<HHH", 0, 1, len(frames)))
    offset = 6 + 16 * len(frames)
    for w, h, data in frames:
        out += struct.pack("<BBBBHHII", w % 256, h % 256, 0, 0, 1, 32, len(data), offset)
        offset += len(data)
    for _, _, data in frames:
        out += data
    path.write_bytes(bytes(out))
    print(f"wrote {path.relative_to(ROOT)} ({len(frames)} sizes, {len(out) // 1024} KB)")


def save_png(im: Image.Image, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    im.save(path, "PNG", optimize=True)
    print(f"wrote {path.relative_to(ROOT)} ({path.stat().st_size // 1024} KB)")


def main() -> None:
    master = load()
    assets, web = ROOT / "assets", ROOT / "web" / "assets"

    save_png(trimmed_logo(master, 600), assets / "logo.png")
    save_png(icon(512, master), assets / "icon.png")
    for s in (16, 24, 32, 48, 64, 128, 256, 512):
        save_png(icon(s, master), assets / "icons" / f"degenprime-{s}.png")

    save_ico([icon(s, master) for s in (16, 20, 24, 32, 40, 48, 64, 96, 128, 256)], assets / "icon.ico")

    icns_sizes = [16, 32, 64, 128, 256, 512, 1024]
    icns = [icon(s, master) for s in icns_sizes]
    icns[-1].save(assets / "icon.icns", format="ICNS", append_images=icns[:-1])
    print("wrote assets/icon.icns")

    logo = trimmed_logo(master, 560)
    logo.save(web / "logo.webp", "WEBP", quality=90, method=6)
    print(f"wrote web/assets/logo.webp ({(web / 'logo.webp').stat().st_size // 1024} KB)")
    save_png(badge(96, master.crop(HEAD), 0.98), web / "mark.png")
    save_ico([icon(s, master) for s in (16, 32, 48)], web / "favicon.ico")
    save_png(icon(192, master), web / "icon-192.png")
    save_png(icon(512, master), web / "icon-512.png")
    save_png(badge(180, master.crop(BUST), 1.0, rounded=False).convert("RGB"), web / "apple-touch-icon.png")


if __name__ == "__main__":
    main()
