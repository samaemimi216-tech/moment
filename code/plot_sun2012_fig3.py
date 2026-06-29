from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DATA_DIRS = [ROOT / "data", ROOT / "三维data", ROOT / "output"]
OUT_DIR = ROOT / "pic"
OUT_FILE = OUT_DIR / "sun2012_fig3_isotherms_600dpi.png"

DPI = 600
FIG_W_PX = 5361
FIG_H_PX = 4474
SINGLE_FIG_W_PX = 4400
SINGLE_FIG_H_PX = 3720

PANELS = [
    ("fig3_t0005.csv", "", "", "8A.png"),
    ("fig3_t0015.csv", "", "", "8B.png"),
    ("fig3_t0050.csv", "", "", "8C.png"),
    ("fig3_steady.csv", "", "", "8D.png"),
]

LEVELS = np.linspace(0.5, 1.0, 21)


def find_data_file(csv_name: str) -> Path:
    for data_dir in DATA_DIRS:
        path = data_dir / csv_name
        if path.exists():
            return path
    searched = ", ".join(str(path / csv_name) for path in DATA_DIRS)
    raise FileNotFoundError(f"Could not find {csv_name}; searched: {searched}")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    names = (
        ("arialbd.ttf", "arial.ttf") if bold else ("arial.ttf", "Arial.ttf")
    )
    for name in names:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default(size=size)


def read_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append((float(row["y"]), float(row["z"]), float(row["theta"])))

    y = np.array(sorted({row[0] for row in rows}), dtype=float)
    z = np.array(sorted({row[1] for row in rows}), dtype=float)
    theta = np.full((z.size, y.size), np.nan, dtype=float)
    y_index = {value: idx for idx, value in enumerate(y)}
    z_index = {value: idx for idx, value in enumerate(z)}
    for yy, zz, tt in rows:
        theta[z_index[zz], y_index[yy]] = tt

    if np.isnan(theta).any():
        raise ValueError(f"Incomplete y-z grid in {path}")
    return y, z, theta


def interpolate_panel(
    y: np.ndarray,
    z: np.ndarray,
    theta: np.ndarray,
    size: int,
) -> np.ndarray:
    xi = np.linspace(y.min(), y.max(), size)
    zi = np.linspace(z.max(), z.min(), size)

    x_idx = np.searchsorted(y, xi, side="right") - 1
    z_idx = np.searchsorted(z, zi, side="right") - 1
    x_idx = np.clip(x_idx, 0, y.size - 2)
    z_idx = np.clip(z_idx, 0, z.size - 2)

    x0 = y[x_idx]
    x1 = y[x_idx + 1]
    z0 = z[z_idx]
    z1 = z[z_idx + 1]
    tx = (xi - x0) / (x1 - x0)
    tz = (zi - z0) / (z1 - z0)

    q00 = theta[z_idx[:, None], x_idx[None, :]]
    q10 = theta[z_idx[:, None], (x_idx + 1)[None, :]]
    q01 = theta[(z_idx + 1)[:, None], x_idx[None, :]]
    q11 = theta[(z_idx + 1)[:, None], (x_idx + 1)[None, :]]

    tx2 = tx[None, :]
    tz2 = tz[:, None]
    values = (
        (1.0 - tx2) * (1.0 - tz2) * q00
        + tx2 * (1.0 - tz2) * q10
        + (1.0 - tx2) * tz2 * q01
        + tx2 * tz2 * q11
    )
    return np.clip(values, LEVELS.min(), LEVELS.max())


def color_ramp(values: np.ndarray) -> np.ndarray:
    # Anchors are sampled from the existing pic/*.png figures so this contour
    # panel sits in the same visual family as the line-plot figures.
    stops = np.array(
        [
            [0.50, 68, 69, 114],
            [0.56, 80, 77, 113],
            [0.62, 127, 200, 210],
            [0.68, 185, 207, 161],
            [0.74, 211, 211, 211],
            [0.80, 204, 173, 204],
            [0.86, 255, 164, 30],
            [0.92, 247, 112, 105],
            [1.00, 197, 85, 94],
        ],
        dtype=float,
    )
    v = np.asarray(values)
    rgb = np.empty((*v.shape, 3), dtype=np.uint8)
    for channel in range(3):
        rgb[..., channel] = np.interp(v, stops[:, 0], stops[:, channel + 1])
    return rgb


def filled_contour_image(values: np.ndarray) -> Image.Image:
    band = np.digitize(values, LEVELS, right=False) - 1
    band = np.clip(band, 0, len(LEVELS) - 2)
    mid = (LEVELS[:-1] + LEVELS[1:]) / 2.0
    colors = color_ramp(mid)[band]
    return Image.fromarray(colors, mode="RGB")


def draw_rotated_text(
    canvas: Image.Image,
    xy: tuple[int, int],
    text: str,
    text_font: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
) -> None:
    bbox = text_font.getbbox(text)
    w = bbox[2] - bbox[0] + 10
    h = bbox[3] - bbox[1] + 10
    layer = Image.new("RGBA", (w, h), (255, 255, 255, 0))
    draw = ImageDraw.Draw(layer)
    draw.text((5 - bbox[0], 5 - bbox[1]), text, font=text_font, fill=fill)
    canvas.alpha_composite(layer.rotate(90, expand=True), xy)


def draw_axes(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    left: int,
    top: int,
    size: int,
    fonts: dict[str, ImageFont.FreeTypeFont],
) -> None:
    ink = (0, 0, 0)
    draw.rectangle([left, top, left + size, top + size], outline=ink, width=10)

    tick_len = 30
    for value in np.linspace(0.0, 1.0, 6):
        x = left + int(round(value * size))
        y = top + size - int(round(value * size))
        label = f"{value:.1f}"

        draw.line([x, top + size, x, top + size + tick_len], fill=ink, width=8)
        box = draw.textbbox((0, 0), label, font=fonts["tick"])
        draw.text(
            (x - (box[2] - box[0]) / 2, top + size + 38),
            label,
            font=fonts["tick"],
            fill=ink,
        )

        draw.line([left - tick_len, y, left, y], fill=ink, width=8)
        box = draw.textbbox((0, 0), label, font=fonts["tick"])
        draw.text(
            (left - tick_len - 30 - (box[2] - box[0]), y - (box[3] - box[1]) / 2),
            label,
            font=fonts["tick"],
            fill=ink,
        )

    y_label = "y*"
    z_label = "z*"
    box = draw.textbbox((0, 0), y_label, font=fonts["label"])
    draw.text(
        (left + size / 2 - (box[2] - box[0]) / 2, top + size + 118),
        y_label,
        font=fonts["label"],
        fill=ink,
    )
    draw_rotated_text(
        canvas,
        (left - 178, top + size // 2 - 48),
        z_label,
        fonts["label"],
        ink,
    )


def draw_colorbar(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    left: int,
    top: int,
    width: int,
    height: int,
    fonts: dict[str, ImageFont.FreeTypeFont],
) -> None:
    values = np.linspace(LEVELS.max(), LEVELS.min(), height)
    bar = Image.fromarray(color_ramp(values).reshape(height, 1, 3), mode="RGB")
    bar = bar.resize((width, height), Image.Resampling.NEAREST)
    canvas.paste(bar, (left, top))

    ink = (0, 0, 0)
    draw.rectangle([left, top, left + width, top + height], outline=ink, width=4)
    for value in np.linspace(0.5, 1.0, 6):
        y = top + height - int(round((value - 0.5) / 0.5 * height))
        draw.line([left + width, y, left + width + 22, y], fill=ink, width=4)
        label = f"{value:.1f}"
        box = draw.textbbox((0, 0), label, font=fonts["tick"])
        draw.text(
            (left + width + 34, y - (box[3] - box[1]) / 2),
            label,
            font=fonts["tick"],
            fill=ink,
        )


def draw_panel(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    csv_name: str,
    panel_label: str,
    time_label: str,
    left: int,
    top: int,
    panel_size: int,
    cbar_width: int,
    cbar_gap: int,
    fonts: dict[str, ImageFont.FreeTypeFont],
) -> None:
    y, z, theta = read_grid(find_data_file(csv_name))
    values = interpolate_panel(y, z, theta, panel_size)
    panel = filled_contour_image(values)
    canvas.paste(panel, (left, top))

    draw_axes(canvas, draw, left, top, panel_size, fonts)
    draw_colorbar(
        canvas,
        draw,
        left + panel_size + cbar_gap,
        top,
        cbar_width,
        panel_size,
        fonts,
    )

    if panel_label:
        panel_x = left - max(108, int(fonts["panel"].size * 1.05))
        panel_y = top + panel_size + max(96, int(fonts["panel"].size * 1.30))
        draw.text(
            (panel_x, panel_y),
            panel_label,
            font=fonts["panel"],
            fill=(35, 35, 35),
        )

    if time_label:
        time_box = draw.textbbox((0, 0), time_label, font=fonts["title"])
        title_y = top - max(98, int(fonts["title"].size * 1.45))
        draw.text(
            (left + panel_size / 2 - (time_box[2] - time_box[0]) / 2, title_y),
            time_label,
            font=fonts["title"],
            fill=(80, 77, 113),
        )


def make_fonts(tick: int, label: int, panel: int, title: int) -> dict[str, ImageFont.FreeTypeFont]:
    return {
        "tick": font(tick),
        "label": font(label, bold=True),
        "panel": font(panel, bold=True),
        "title": font(title, bold=True),
    }


def write_composite() -> Path:
    canvas = Image.new("RGBA", (FIG_W_PX, FIG_H_PX), "white")
    draw = ImageDraw.Draw(canvas)

    fonts = make_fonts(tick=64, label=86, panel=72, title=70)

    panel_size = 1320
    cbar_width = 82
    cbar_gap = 66
    lefts = [420, 2925]
    tops = [350, 2325]

    for index, (csv_name, panel_label, time_label, _out_name) in enumerate(PANELS):
        row, col = divmod(index, 2)
        left = lefts[col]
        top = tops[row]
        draw_panel(
            canvas,
            draw,
            csv_name,
            panel_label,
            time_label,
            left,
            top,
            panel_size,
            cbar_width,
            cbar_gap,
            fonts,
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rgb = canvas.convert("RGB")
    rgb.save(OUT_FILE, dpi=(DPI, DPI), optimize=True)
    return OUT_FILE


def write_single_panels() -> list[Path]:
    outputs: list[Path] = []
    fonts = make_fonts(tick=118, label=158, panel=132, title=128)
    panel_size = 2920
    cbar_width = 150
    cbar_gap = 120
    left = 620
    top = 230

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for csv_name, panel_label, time_label, out_name in PANELS:
        canvas = Image.new("RGBA", (SINGLE_FIG_W_PX, SINGLE_FIG_H_PX), "white")
        draw = ImageDraw.Draw(canvas)
        draw_panel(
            canvas,
            draw,
            csv_name,
            panel_label,
            time_label,
            left,
            top,
            panel_size,
            cbar_width,
            cbar_gap,
            fonts,
        )
        output = OUT_DIR / out_name
        canvas.convert("RGB").save(output, dpi=(DPI, DPI), optimize=True)
        outputs.append(output)
    return outputs


def main() -> None:
    outputs = [write_composite(), *write_single_panels()]
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
