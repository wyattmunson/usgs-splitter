#!/usr/bin/env python3
"""
usgs_splitter.py

Split a USGS quadrangle GeoTIFF into four 8.5×11 printable pages
(NW, NE, SW, SE) plus a cover page.

  - Auto-detects and crops the neatline (removes title/legend margins)
  - Cover page includes the three collar elements from the original sheet:
    state locator mini-map, adjoining quadrangles index, GN/MN declination

Usage:
    python usgs_splitter.py input.tif [-o output.pdf] [--name "Quad Name"]
"""

import sys
import os
import math
import argparse
import re

import numpy as np
import rasterio
from rasterio.warp import transform_bounds
from rasterio.windows import Window
from rasterio.enums import Resampling
from pyproj import Transformer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def dms_str(degrees, is_lat):
    """Convert decimal degrees to a DMS string like 37°30'15.0\"N."""
    neg = degrees < 0
    d   = int(abs(degrees))
    rem = (abs(degrees) - d) * 60
    m   = int(rem)
    s   = (rem - m) * 60
    suffix = ('S' if neg else 'N') if is_lat else ('W' if neg else 'E')
    return f"{d}°{m:02d}'{s:04.1f}\"{suffix}"


def nice_ticks(lo, hi, target_n=5):
    """Return evenly-spaced 'nice' tick values spanning [lo, hi]."""
    span = hi - lo
    candidates = [
        1/720, 1/360, 1/180, 1/120, 1/60, 1/30,
        1/20,  1/12,  1/8,   1/6,   1/4,  1/3,  1/2, 1.0,
    ]
    step  = min(candidates, key=lambda s: abs(span / s - target_n))
    start = math.ceil(lo / step) * step
    ticks, v = [], start
    while v <= hi + 1e-9:
        ticks.append(round(v, 9))
        v += step
    return ticks


# ---------------------------------------------------------------------------
# Neatline detection
# ---------------------------------------------------------------------------

def find_neatline(src):
    """
    Locate the map neatline in pixel coordinates by snapping the file's
    geographic extent inward to the nearest USGS 7.5-minute grid (0.125°),
    then projecting those four corners back to pixel space.

    Returns a dict with:
      top, bottom, left, right  – full-resolution pixel bounds of the map area
      quad_bounds               – (west, south, east, north) WGS84
    """
    crs = src.crs or "EPSG:4326"
    file_w, file_s, file_e, file_n = transform_bounds(crs, "EPSG:4326", *src.bounds)

    # Snap inward to the nearest 7.5-minute boundary
    step   = 0.125
    quad_n = math.floor(file_n / step) * step
    quad_s = math.ceil (file_s / step) * step
    quad_w = math.ceil (file_w / step) * step   # ceil of negative → less negative (east)
    quad_e = math.floor(file_e / step) * step   # floor of negative → more negative (west)

    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    corners_native = [
        tr.transform(lon, lat)
        for lon, lat in [
            (quad_w, quad_n), (quad_e, quad_n),
            (quad_w, quad_s), (quad_e, quad_s),
        ]
    ]
    inv    = ~src.transform
    pixels = [inv * pt for pt in corners_native]
    rows   = [r for c, r in pixels]
    cols   = [c for c, r in pixels]

    return {
        'top':         max(0, round(min(rows))),
        'bottom':      min(src.height, round(max(rows))),
        'left':        max(0, round(min(cols))),
        'right':       min(src.width,  round(max(cols))),
        'quad_bounds': (quad_w, quad_s, quad_e, quad_n),
    }


# ---------------------------------------------------------------------------
# Collar element extraction
# ---------------------------------------------------------------------------

# Horizontal positions of the three elements expressed as fractions of image
# width.  These match the USGS US Topo standard bottom-collar layout:
#
#   0%–27%   production credits (text block)
#   27%–42%  declination diagram  (GN / MN arrows)
#   42%–59%  scale bar
#   59%–74%  state locator + adjoining quadrangles
#   74%–84%  road classification
#   84%–100% title / barcode
#
# The top ~100 px of the collar are coordinate tick marks; we skip them.

_COLLAR_ELEMENTS = {
    # key: (x_start_frac, x_width_frac, y_rel_start_frac, y_rel_height_frac)
    # y_rel fractions are of the usable collar height (after the 100-px skip)
    #
    # Horizontal layout of the USGS US Topo bottom collar:
    #   0–27%   production credits
    #   27–42%  declination diagram (GN/MN arrows + labels)   ← isolated here
    #   42–59%  scale bar
    #   59–74%  state locator + adjoining quadrangles          ← isolated here
    #   74–84%  road classification
    #   84–100% title / barcode
    'declination':   (0.305, 0.078, 0.00, 0.55),  # stop before the US National Grid box
    'state_locator': (0.620, 0.090, 0.00, 0.50),
    'adj_quads':     (0.595, 0.160, 0.50, 0.38),  # tighter — grid+names only
}


def _extract_collar_elements(src, neatline):
    """
    Read the three useful diagrams from the bottom collar of the USGS sheet.
    Returns a dict of uint8 ndarray (H, W, 3) or None if not extractable.
    """
    H, W      = src.height, src.width
    collar_y0 = neatline['bottom'] + 100   # skip coordinate tick-mark row
    collar_h  = H - collar_y0

    if collar_h < 50:
        return {k: None for k in _COLLAR_ELEMENTS}

    result = {}
    for key, (xf, wf, yf, hf) in _COLLAR_ELEMENTS.items():
        x0  = int(W * xf)
        ew  = int(W * wf)
        y0  = collar_y0 + int(collar_h * yf)
        eh  = int(collar_h * hf)
        eh  = min(eh, H - y0)
        if eh <= 0 or ew <= 0:
            result[key] = None
            continue
        win  = Window(x0, y0, ew, eh)
        data = src.read(window=win)
        if data.shape[0] == 1:
            rgb = np.stack([data[0]] * 3, axis=-1)
        elif data.shape[0] >= 3:
            rgb = np.transpose(data[:3], (1, 2, 0))
        else:
            rgb = np.transpose(np.concatenate([data, data[:1]], axis=0)[:3], (1, 2, 0))
        result[key] = rgb.astype(np.uint8)

    return result


def _autocrop_white(arr, pad=12):
    """Trim excess white border from a uint8 HxWx3 array."""
    if arr is None:
        return None
    diff = np.abs(arr.astype(int) - 255).sum(axis=2)
    rows = np.where(diff.max(axis=1) > 8)[0]
    cols = np.where(diff.max(axis=0) > 8)[0]
    if len(rows) == 0 or len(cols) == 0:
        return arr
    r0 = max(0, rows[0]  - pad)
    r1 = min(arr.shape[0], rows[-1] + pad + 1)
    c0 = max(0, cols[0]  - pad)
    c1 = min(arr.shape[1], cols[-1] + pad + 1)
    return arr[r0:r1, c0:c1]


# ---------------------------------------------------------------------------
# Image loading
# ---------------------------------------------------------------------------

def load_geotiff(path, max_pixels=30_000_000):
    """
    Open a USGS quadrangle GeoTIFF, crop it to the map neatline, and
    extract the three collar diagrams.

    Returns:
      img         – uint8 ndarray (H, W, 3) cropped to the map area
      bounds      – (west, south, east, north) WGS84 of the cropped region
      collar_imgs – dict: 'declination', 'state_locator', 'adj_quads'
    """
    with rasterio.open(path) as src:
        neatline = find_neatline(src)
        top, bottom, left, right = (
            neatline['top'], neatline['bottom'],
            neatline['left'], neatline['right'],
        )
        map_h = bottom - top
        map_w = right  - left

        scale = 1.0
        if map_h * map_w > max_pixels:
            scale = math.sqrt(max_pixels / (map_h * map_w))
        out_h = max(1, int(map_h * scale))
        out_w = max(1, int(map_w * scale))

        win  = Window(left, top, map_w, map_h)
        data = src.read(
            window=win,
            out_shape=(src.count, out_h, out_w),
            resampling=Resampling.lanczos,
        )
        collar_imgs = _extract_collar_elements(src, neatline)

    # Build RGB
    bands = data.shape[0]
    if bands == 1:
        rgb = np.stack([data[0]] * 3, axis=-1)
    elif bands >= 3:
        rgb = np.transpose(data[:3], (1, 2, 0))
    else:
        padded = np.concatenate([data, data[:1]], axis=0)
        rgb    = np.transpose(padded[:3], (1, 2, 0))

    if rgb.dtype != np.uint8:
        lo, hi_v = float(rgb.min()), float(rgb.max())
        if hi_v > lo:
            rgb = ((rgb.astype(np.float32) - lo) / (hi_v - lo) * 255).astype(np.uint8)
        else:
            rgb = np.zeros_like(rgb, dtype=np.uint8)

    # Auto-crop any white fringe left by the neatline snap
    collar_imgs = {k: _autocrop_white(v) for k, v in collar_imgs.items()}

    return rgb, neatline['quad_bounds'], collar_imgs


# ---------------------------------------------------------------------------
# Quadrant splitting
# ---------------------------------------------------------------------------

QUAD_ORDER  = ['NW', 'NE', 'SW', 'SE']
QUAD_PAGES  = {'NW': 2, 'NE': 3, 'SW': 4, 'SE': 5}
QUAD_COLORS = {'NW': '#d4e6f1', 'NE': '#d5f5e3', 'SW': '#fdebd0', 'SE': '#e8daef'}


def split_quadrants(img, bounds):
    """Split the cropped map array into NW / NE / SW / SE quadrants."""
    h, w   = img.shape[:2]
    mh, mw = h // 2, w // 2
    west, south, east, north = bounds
    mid_lat = (north + south) / 2
    mid_lon = (east  + west)  / 2

    return {
        'NW': (img[:mh, :mw],  (west,    mid_lat, mid_lon, north)),
        'NE': (img[:mh, mw:],  (mid_lon, mid_lat, east,    north)),
        'SW': (img[mh:, :mw],  (west,    south,   mid_lon, mid_lat)),
        'SE': (img[mh:, mw:],  (mid_lon, south,   east,    mid_lat)),
    }


# ---------------------------------------------------------------------------
# Quadrant page renderer
# ---------------------------------------------------------------------------

def render_quadrant_page(pdf, name, img, bounds, page_num, map_name, dpi):
    """
    Render one quadrant as a portrait 8.5 × 11 inch page.

    The image is displayed in its native pixel space (no extent= argument),
    so matplotlib's default square-pixel aspect is used and nothing is
    stretched.  Geographic coordinate labels are placed by converting lon/lat
    values to pixel positions manually.
    """
    west, south, east, north = bounds
    h_px, w_px = img.shape[:2]
    mid_lat = (north + south) / 2

    fig, ax = plt.subplots(figsize=(8.5, 11), dpi=dpi)
    fig.patch.set_facecolor('white')

    # Show in pixel space — imshow default is aspect='equal' (square pixels),
    # which preserves the image's true proportions.
    ax.imshow(img, origin='upper', interpolation='bilinear')
    # imshow sets xlim=(-0.5, w_px-0.5), ylim=(h_px-0.5, -0.5)

    # ── Geographic ticks mapped to pixel positions ─────────────────────
    def col(lon): return (lon - west) / (east - west) * w_px - 0.5
    def row(lat): return (north - lat) / (north - south) * h_px - 0.5

    lon_ticks = nice_ticks(west, east, target_n=5)
    lat_ticks = nice_ticks(south, north, target_n=6)

    ax.set_xticks([col(l) for l in lon_ticks])
    ax.set_xticklabels([dms_str(l, False) for l in lon_ticks])
    ax.set_yticks([row(l) for l in lat_ticks])
    ax.set_yticklabels([dms_str(l, True) for l in lat_ticks])

    ax.tick_params(
        axis='both', which='major',
        top=True, bottom=True, left=True, right=True,
        labeltop=True, labelbottom=True, labelleft=True, labelright=True,
        labelsize=6, length=5, width=0.8,
    )
    ax.tick_params(axis='x', labelrotation=35)
    ax.grid(True, color='#aaaaaa', linestyle='--', linewidth=0.35, alpha=0.55)

    # ── Scale bar (pixel coordinates) ─────────────────────────────────
    miles_per_deg_lon = 69.0 * math.cos(math.radians(mid_lat))
    target_miles = (east - west) * miles_per_deg_lon * 0.20
    nice_m    = [0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 25, 50, 100]
    bar_miles = min(nice_m, key=lambda x: abs(x - target_miles))
    bar_km    = bar_miles * 1.60934
    bar_px    = bar_miles / miles_per_deg_lon / (east - west) * w_px

    sb_x0 = w_px * 0.04
    sb_y0 = h_px * 0.968
    sb_x1 = sb_x0 + bar_px
    sb_th = h_px * 0.007

    ax.plot([sb_x0, sb_x1], [sb_y0, sb_y0], color='black',
            linewidth=3, solid_capstyle='butt', zorder=5)
    for x in (sb_x0, sb_x1):
        ax.plot([x, x], [sb_y0 - sb_th, sb_y0 + sb_th],
                color='black', linewidth=1.5, zorder=5)
    ax.text((sb_x0 + sb_x1) / 2, sb_y0 - sb_th * 1.5,
            f'{bar_miles:g} mi  /  {bar_km:.1f} km',
            ha='center', va='bottom', fontsize=6.5, zorder=6,
            bbox=dict(facecolor='white', alpha=0.75, edgecolor='none', pad=1))

    # ── Quadrant badge ─────────────────────────────────────────────────
    ax.text(w_px * 0.012, h_px * 0.020,
            f'{name}  —  {map_name}',
            ha='left', va='top', fontsize=8, fontweight='bold', color='white',
            zorder=7,
            bbox=dict(facecolor='#0050a0', alpha=0.80,
                      edgecolor='none', boxstyle='round,pad=0.35'))

    fig.text(0.985, 0.012, f'Page {page_num}',
             ha='right', va='bottom', fontsize=8, color='white', bbox=dict(facecolor='#a00000', alpha=0.80,
                      edgecolor='none', boxstyle='round,pad=0.35'))

    plt.tight_layout(pad=0.55)
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Cover page renderer
# ---------------------------------------------------------------------------

def render_cover_page(pdf, bounds, map_name, quadrants, collar_imgs):
    """
    Render the cover page (portrait 8.5 × 11).

    Layout (top to bottom):
      Title block
      ── blue rule ──
      [Declination] | [State locator] | [Adjoining quads]   ← collar panels
      ── 2×2 layout diagram ──
      Print instructions
      Footer
    """
    west, south, east, north = bounds
    mid_lat = (north + south) / 2
    mid_lon = (east  + west)  / 2

    fig = plt.figure(figsize=(8.5, 11), dpi=300)
    fig.patch.set_facecolor('white')

    # ── Title ────────────────────────────────────────────────────────────
    fig.text(0.5, 0.970, 'USGS Quadrangle  ·  Print Split',
             ha='center', va='top', fontsize=12, fontweight='bold', color='#222')
    fig.text(0.5, 0.934, map_name,
             ha='center', va='top', fontsize=16, color='#0050a0', fontweight='bold')
    fig.text(0.5, 0.902,
             f'{dms_str(west, False)} – {dms_str(east, False)}'
             f'   ·   {dms_str(south, True)} – {dms_str(north, True)}',
             ha='center', va='top', fontsize=8, color='#555')

    rule_ax = fig.add_axes([0.08, 0.891, 0.84, 0.003])
    rule_ax.axhline(0.5, color='#0050a0', linewidth=1.5)
    rule_ax.axis('off')

    # ── Collar panels ────────────────────────────────────────────────────
    # Each entry: (collar_key, display_label, panel_width_frac, panel_height_frac)
    # adj_quads gets a wider box to suit its landscape 3×3 grid layout.
    PANEL_GAP = 0.018
    panel_defs = [
        ('declination',   'Grid & Magnetic North', 0.240, 0.185),
        ('state_locator', 'Quadrangle Location',   0.240, 0.185),
        ('adj_quads',     'Adjoining Quadrangles', 0.310, 0.185),
    ]
    total_w    = sum(w for _, _, w, _ in panel_defs) + PANEL_GAP * (len(panel_defs) - 1)
    panel_left = (1.0 - total_w) / 2
    panel_top  = 0.881   # top edge of panels in figure coords

    x_cursor = panel_left
    for key, label, pw, ph in panel_defs:
        y0   = panel_top - ph
        ax_p = fig.add_axes([x_cursor, y0, pw, ph])
        ax_p.set_facecolor('#f9f9f9')
        for sp in ax_p.spines.values():
            sp.set_edgecolor('#cccccc')
            sp.set_linewidth(0.8)
        ax_p.set_xticks([])
        ax_p.set_yticks([])

        arr = collar_imgs.get(key)
        if arr is not None and arr.size > 0:
            # aspect='auto' fills the panel box — avoids blank space when the
            # image's pixel ratio doesn't match the panel shape.
            ax_p.imshow(arr, aspect='auto', interpolation='bilinear')
        else:
            ax_p.text(0.5, 0.5, '—', ha='center', va='center',
                      fontsize=10, color='#bbb', transform=ax_p.transAxes)

        ax_p.set_title(label, fontsize=7, color='#555', pad=3)
        x_cursor += pw + PANEL_GAP

    # ── 2×2 layout diagram ───────────────────────────────────────────────
    max_panel_h = max(ph for _, _, _, ph in panel_defs)
    diag_top    = panel_top - max_panel_h - 0.022
    diag_height = 0.365
    ax = fig.add_axes([0.10, diag_top - diag_height, 0.80, diag_height])
    ax.set_xlim(-0.18, 2.18)
    ax.set_ylim(-0.18, 2.18)
    ax.set_aspect('equal')
    ax.axis('off')

    grid = [('NW', 0, 1), ('NE', 1, 1), ('SW', 0, 0), ('SE', 1, 0)]
    for qname, col, row in grid:
        page  = QUAD_PAGES[qname]
        color = QUAD_COLORS[qname]
        rect  = mpatches.FancyBboxPatch(
            (col + 0.04, row + 0.04), 0.92, 0.92,
            boxstyle='round,pad=0.02',
            facecolor=color, edgecolor='#334', linewidth=1.5,
        )
        ax.add_patch(rect)
        ax.text(col + 0.50, row + 0.64, qname,
                ha='center', va='center', fontsize=22, fontweight='bold', color='#222')
        ax.text(col + 0.50, row + 0.30, f'Page {page}',
                ha='center', va='center', fontsize=11, color='#444')

    for txt, x, y in [('N', 1.0, 2.14), ('S', 1.0, -0.14),
                       ('W', -0.14, 1.0), ('E', 2.14, 1.0)]:
        ax.text(x, y, txt, ha='center', va='center',
                fontsize=12, fontweight='bold', color='#444')

    ax.axhline(1, xmin=0.02, xmax=0.98, color='#777', linewidth=0.9, linestyle='--')
    ax.axvline(1, ymin=0.02, ymax=0.98, color='#777', linewidth=0.9, linestyle='--')
    ax.text(1.0,   1.048, dms_str(mid_lat, True),
            ha='center', va='bottom', fontsize=7, color='#555')
    ax.text(1.048, 1.0,   dms_str(mid_lon, False),
            ha='left', va='center', fontsize=7, color='#555', rotation=90)

    fs = 6.5
    for x, y, ha, va, lat, lon in [
        (0.02, 1.98, 'left',  'top',    north, west),
        (1.98, 1.98, 'right', 'top',    north, east),
        (0.02, 0.02, 'left',  'bottom', south, west),
        (1.98, 0.02, 'right', 'bottom', south, east),
    ]:
        ax.text(x, y, f'{dms_str(lat, True)}\n{dms_str(lon, False)}',
                ha=ha, va=va, fontsize=fs, color='#555', linespacing=1.4)

    # ── Instructions ─────────────────────────────────────────────────────
    instr_top = diag_top - diag_height - 0.016
    fig.text(0.5, instr_top, 'Print Instructions',
             ha='center', va='top', fontsize=10, fontweight='bold', color='#333')

    lines = [
        ('Page 1  (this cover)',        'portrait orientation'),
        ('Pages 2–5  (quadrant maps)',  'landscape orientation, 100 % scale'),
        ('Arrange printed pages',       'match the layout diagram above'),
        ('Adjacent pages share',        'a common map edge; coordinate labels align'),
    ]
    for i, (label, detail) in enumerate(lines):
        y = instr_top - 0.028 - i * 0.028
        fig.text(0.14, y, f'•  {label}',
                 ha='left', va='top', fontsize=8.5, color='#222', fontweight='bold')
        fig.text(0.50, y, detail,
                 ha='left', va='top', fontsize=8.5, color='#555')

    tip_y = instr_top - 0.028 - 4 * 0.028
    fig.text(0.14, tip_y,
             'Tip:  Choose "Actual Size" in your print dialog to preserve scale bar accuracy.',
             ha='left', va='top', fontsize=8, color='#555', style='italic')

    # ── Footer ───────────────────────────────────────────────────────────
    fig.text(0.5, 0.022,
             'As with all maps, inaccuracies may exist and conditions may change. User assumes all risk associated with the use of this map.\nGenerated by usgs_splitter.py  ·  USGS data is in the public domain.',
             ha='center', va='bottom', fontsize=7, color='#999', style='italic')

    pdf.savefig(fig)
    plt.close(fig)



# ---------------------------------------------------------------------------
# Name formatter
# ---------------------------------------------------------------------------
def derive_map_name(input_path):
    """
    Returns a formatted map name derived from the file name. Returns a 
    name like `Sea Cliff, NY`. This assumes a city name does not contain
    a number and states are always two charachters. 
    
    Args:
        input_path (str): File path of .tif file

    Returns:
        str: Formatted name 
    """
    base = os.path.splitext(os.path.basename(input_path))[0]
    parts = base.split('_')
    
    state = parts[0]
    # find index of first part that looks like a number (the date)
    date_idx = next((i for i, p in enumerate(parts) if p.isdigit()), len(parts))
    
    city_parts = parts[1:date_idx]
    city = ' '.join(p.title() for p in city_parts)
    
    return f"{city}, {state}"

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Split a USGS quadrangle GeoTIFF into 4 printable 8.5×11 pages.'
    )
    parser.add_argument('input',
                        help='Path to the GeoTIFF file')
    parser.add_argument('-o', '--output',
                        help='Output PDF path (default: outputs/<input>_split.pdf)')
    parser.add_argument('--name',
                        help='Quadrangle name for display (default: from filename)')
    parser.add_argument('--max-pixels', type=int, default=30_000_000, metavar='N',
                        help='Downsample if the map area exceeds N pixels (default: 30 000 000)')
    parser.add_argument('--dpi', '-d',
                        help='DPI of split maps (default: 300 DPI)')
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f'Error: file not found: {args.input}')

    if args.output:
        output_path = args.output
    else:
        repo_root  = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(repo_root, 'outputs')
        os.makedirs(output_dir, exist_ok=True)
        base_name  = os.path.splitext(os.path.basename(args.input))[0]
        output_path = os.path.join(output_dir, base_name + '_split.pdf')
    # map_name    = (args.name
    #                or os.path.splitext(os.path.basename(args.input))[0]
    #                   .replace('_', ' ').replace('-', ' ').title())
    map_name = args.name or derive_map_name(args.input)
    dpi = int(args.dpi or 300)

    print(f'Loading   {args.input} …')
    img, bounds, collar_imgs = load_geotiff(args.input, max_pixels=args.max_pixels)
    west, south, east, north = bounds
    print(f'  Extent  : {dms_str(west, False)} → {dms_str(east, False)}'
          f',  {dms_str(south, True)} → {dms_str(north, True)}')
    print(f'  Map size: {img.shape[1]} × {img.shape[0]} px  (neatline-cropped)')
    found = [k for k, v in collar_imgs.items() if v is not None]
    print(f'  Collar elements found: {found}')

    print('Splitting into quadrants …')
    quadrants = split_quadrants(img, bounds)

    print(f'Writing   {output_path} …')
    with PdfPages(output_path) as pdf:
        meta = pdf.infodict()
        meta['Title']   = f'USGS Quadrangle Split: {map_name}'
        meta['Author']  = 'usgs_splitter'
        meta['Subject'] = 'USGS printable quadrangle split'

        render_cover_page(pdf, bounds, map_name, quadrants, collar_imgs)
        print('  [1/5] Cover page')

        for i, name in enumerate(QUAD_ORDER, start=2):
            img_q, bounds_q = quadrants[name]
            render_quadrant_page(pdf, name, img_q, bounds_q, QUAD_PAGES[name], map_name, dpi)
            print(f'  [{i}/5] {name} quadrant  (page {QUAD_PAGES[name]})')

    print(f'\nDone → {output_path}')


if __name__ == '__main__':
    main()
