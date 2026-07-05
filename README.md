# USGS Splitter

A command-line tool that takes a USGS 7.5-minute quadrangle GeoTIFF and splits it into a printable 5-page PDF — a cover page plus four map quadrants (NW, NE, SW, SE) sized for standard 8.5 × 11" letter paper.

- [How it works](#how-it-works)
- [Downloading USGS Quadrangles](#downloading-usgs-quadrangles)
- [Installation](#installation)
- [Usage](#usage)
- [Printing](#printing)

---

## How it works

USGS US Topo quadrangles are large, single-sheet maps. This tool:

1. **Auto-detects the neatline** — snaps the file's geographic bounds to the nearest 7.5-minute grid (multiples of 0.125°) and uses the embedded GeoTIFF coordinate transform to find the exact pixel boundary of the map area, trimming the title block, legend collar, and white margins automatically.
2. **Splits the map** into four equal quadrants at the geographic center point (NW / NE / SW / SE).
3. **Renders each quadrant** to a portrait 8.5 × 11" page at the image's natural aspect ratio (no stretching), with DMS coordinate labels on all four edges, a dashed lat/lon grid, and a scale bar.
4. **Generates a cover page** (portrait) showing:
   - The quadrangle name and full extent in DMS coordinates
   - Three collar diagrams extracted from the original sheet: the GN/MN declination diagram, the state locator mini-map, and the adjoining quadrangles index
   - A 2 × 2 layout diagram with page numbers and corner coordinates
   - Print instructions

### Output PDF structure

| Page | Content                                                     | Orientation |
| ---- | ----------------------------------------------------------- | ----------- |
| 1    | Cover — layout diagram, collar elements, print instructions | Portrait    |
| 2    | NW quadrant                                                 | Portrait    |
| 3    | NE quadrant                                                 | Portrait    |
| 4    | SW quadrant                                                 | Portrait    |
| 5    | SE quadrant                                                 | Portrait    |

---

## Downloading USGS quadrangles

1. Go to the [USGS National Map Downloader](https://apps.nationalmap.gov/downloader/)
   - Or use [USGS topoView](https://ngmdb.usgs.gov/topoview/viewer)
2. Search for a location and select **US Topo** as the product type
3. Download the ZIP — it contains several files including two `.tif` files:
   - `*_TM_geo.tif` — the standard topographic map (use this one)
   - `*_TMorth_geo.tif` — an orthophoto-backed variant (also works)
4. The `.tfw`, `.tif.prj`, and `.tif.xml` sidecar files are **not needed** — all coordinate metadata is embedded in the GeoTIFF itself

> [!NOTE]
> This tool is only tested for 2022 USGS quadrangles (1:24,000). Using prior years or other sizes may not render properly.

---

## Installation

**Python 3.9 or newer is required.**

### Option A — conda (recommended, handles GDAL automatically)

```bash
conda install -c conda-forge rasterio numpy matplotlib pyproj
```

### Option B — pip (requires GDAL to be installed on the system)

```bash
# macOS (Homebrew)
brew install gdal
pip install rasterio numpy matplotlib pyproj

# Ubuntu / Debian
sudo apt install gdal-bin libgdal-dev
pip install rasterio numpy matplotlib pyproj
```

Install from `requirements.txt`:

```bash
pip install -r requirements.txt
```

---

## Usage

```
python usgs_splitter.py <input.tif> [-o <output.pdf>] [--name "Quad Name"] [--max-pixels N]
```

### Arguments

| Argument         | Description                                                                                                                                 |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `input`          | Path to the USGS quadrangle GeoTIFF file **(required)**                                                                                     |
| `-o`, `--output` | Output PDF path. Defaults to `outputs/<input>_split.pdf` in the repo root                                                                  |
| `--name`         | Display name for the quadrangle. Defaults to a cleaned-up version of the filename                                                           |
| `--max-pixels N` | Downsample the map area if it exceeds N pixels before splitting. Useful for limiting memory use on large files. Default: `30000000` (30 MP) |
| `--dpi N`, `-d N` | Split map DPI as integer. Lower DPIs will increase processing time. Input DPI is likely no more than 300. Default: `300`. |

### Examples

Basic usage — output PDF written to `outputs/VA_Richmond_20220920_TM_geo_split.pdf`:

```bash
python usgs_splitter.py VA_Richmond_20220920_TM_geo.tif
```

Specify the output path and display name:

```bash
python usgs_splitter.py VA_Richmond_20220920_TM_geo.tif \
  -o ~/Desktop/Richmond_split.pdf \
  --name "Richmond, Virginia"
```

Limit memory use for a large file:

```bash
python usgs_splitter.py large_quad.tif --max-pixels 15000000
```

### Bulk Processor

Automatically generate multiple maps. Add maps to the `inputs/` directory. Map names will be automatically generated and outputted to the `outputs/` directory. Inputs that have been processed will be moved to the `inputs-completed/` directory.

```bash
./run_batch.sh
```

---

## Printing

- **Page 1 (cover):** print portrait
- **Pages 2–5 (quadrant maps):** print portrait at **100% / Actual Size** — do not use "Fit to Page", as this will shrink the map and throw off the scale bar
- Arrange the four printed pages in the 2 × 2 layout shown on the cover; adjacent pages share a common map edge and coordinate labels will align

---

## Dependencies

| Package      | Purpose                                                         |
| ------------ | --------------------------------------------------------------- |
| `rasterio`   | Reads GeoTIFF files and handles coordinate transforms           |
| `numpy`      | Array manipulation for image data                               |
| `matplotlib` | Renders map pages and generates the PDF                         |
| `pyproj`     | Converts between projected CRS and WGS84 for neatline detection |

All are available on conda-forge and PyPI. See `requirements.txt` for pinned minimum versions.
