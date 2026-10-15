"""Portable font selection for the Pillow structure and orbital figures."""
import os
from PIL import ImageFont


def font(size, bold=False):
    configured = os.environ.get("QUANTUM2_FONT_BOLD" if bold else "QUANTUM2_FONT")
    names = [configured] if configured else [
        "/System/Library/Fonts/Supplemental/Arial" + (" Bold" if bold else "") + ".ttf",
        "C:/Windows/Fonts/" + ("arialbd.ttf" if bold else "arial.ttf"),
        "LiberationSans-Bold.ttf" if bold else "LiberationSans-Regular.ttf",
        "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
    ]
    for name in names:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            pass
    raise RuntimeError("No TrueType font found; set QUANTUM2_FONT and QUANTUM2_FONT_BOLD to font files.")
