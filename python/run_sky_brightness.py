"""
run_sky_brightness.py  —  Distance Project

Batch-processes CR2 fisheye images through DiCaLum and outputs:
  data/sky_brightness_per_image.csv       — one row per image
  data/sky_brightness_measurements.csv    — one row per site (p1/p2/p3 averaged),
     ready to join into the R pipeline via site column.

Run from the Distance_Project root:
    cd /path/to/Distance_Project
    source ../night_sky_project/environment/bin/activate
    python python/run_sky_brightness.py

── CAMERA / LENS ─────────────────────────────────────────────────────────────
cameras: 0=EOS6D  1=EOS60D  2=M100  3=M200  4=SonyA7S
lenses:  0=Sigma8mm  1=Sigma4.5mm  2=Samyang8mm  3=Samyang24mm
         4=Samyang50mm  5=Meike6.5mm  6=NoCorrection
──────────────────────────────────────────────────────────────────────────────
"""

CAMERA_IDX = 0   # EOS6D   ← change if different from ColterBay
LENS_IDX   = 0   # Sigma8mm ← change if different from ColterBay

# ── FILENAME CONVENTION ───────────────────────────────────────────────────────
# Expected pattern: grte_<SITE>_dark_p<N>.CR2
# e.g. grte_grte01_dark_p1.CR2  →  site = "GRTE01"
# If your filenames differ, update extract_site() below.
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
import pandas as pd
from pathlib import Path

_PROJECT_ROOT = Path(__file__).parent.parent
RAW_DIR = Path(
    "/Users/nanderson/Library/CloudStorage/"
    "GoogleDrive-nigel_anderson@brown.edu/"
    ".shortcut-targets-by-id/"
    "1sSdpOAdUOgAVbJGpTKgB3-CJAfnjsvKJ/"
    "grandteton_distanceproject/"
    "grte_distance_nightsky_pics"
)
OUT_DIR = _PROJECT_ROOT / "data"

try:
    import dicalum
except ImportError:
    raise ImportError(
        "DiCaLum not installed.\n"
        "Activate the venv and run:  python -m pip install dicalum==4.0b7"
    )

try:
    import exifread
except ImportError:
    raise ImportError("exifread not installed.  Run:  python -m pip install exifread")


def read_exif(filepath: Path):
    """Read ISO, aperture, shutter from EXIF. Returns (iso, aperture, shutter) as floats."""
    with open(filepath, "rb") as f:
        tags = exifread.process_file(f)

    def parse(tag):
        val = str(tags.get(tag, "0"))
        if "/" in val:
            n, d = val.split("/")
            return float(n) / float(d)
        return float(val) if val else 0.0

    return parse("EXIF ISOSpeedRatings"), parse("EXIF FNumber"), parse("EXIF ExposureTime")


# Hide the Tk window DiCaLum opens at import
dicalum.TopWin.withdraw()

dicalum.dclinst.camera = CAMERA_IDX
dicalum.dclinst.lens   = LENS_IDX
dicalum.setvig()
print(f"Camera : {dicalum.cameras[CAMERA_IDX]}")
print(f"Lens   : {dicalum.lenses[LENS_IDX]}")
print(f"Vignetting matrix set — {dicalum.dcldat.V.shape}\n")


def parse_filename(filename_stem: str) -> tuple[str, str]:
    """
    'grte01_distance_dark_p1'    → ('GRTE01', 'dark')
    'grte01_distance_white100_p1' → ('GRTE01', 'white100')
    Returns (site, condition).
    """
    parts = filename_stem.split("_")
    site      = parts[0].upper() if parts else filename_stem.upper()
    condition = parts[2].lower() if len(parts) >= 3 else "unknown"
    return site, condition


def process_file(cr2_path: Path) -> dict | None:
    print(f"  Processing: {cr2_path.name} ...", end=" ", flush=True)

    iso, aperture, shutter = read_exif(cr2_path)

    if aperture == 0:
        aperture = dicalum.LensList[LENS_IDX].aper

    if iso == 0 or shutter == 0:
        print(f"FAILED (EXIF missing — iso={iso} t={shutter})")
        return None

    ddat = dicalum.rawread(str(cr2_path))

    if isinstance(ddat, int) and ddat == -1:
        print("FAILED (rawread returned -1)")
        return None

    try:
        epo     = 6400 * aperture ** 2 / (iso * shutter)
        r_dsu   = epo * ddat.R
        g_dsu   = epo * ddat.G
        b_dsu   = epo * ddat.B
        lum_dsu = 0.2126 * r_dsu + 0.7152 * g_dsu + 0.0722 * b_dsu

        sky = ddat.M

        mean_dsu   = float(np.mean(lum_dsu[sky]))
        median_dsu = float(np.median(lum_dsu[sky]))
        sd_dsu     = float(np.std(lum_dsu[sky]))
        mean_dsu_g = float(np.mean(g_dsu[sky]))

        print(f"done  (mean DSU = {mean_dsu:.4f})")

        site, condition = parse_filename(cr2_path.stem)

        return {
            "filename":   cr2_path.name,
            "site":       site,
            "condition":  condition,
            "mean_dsu":   mean_dsu,
            "median_dsu": median_dsu,
            "sd_dsu":     sd_dsu,
            "mean_dsu_g": mean_dsu_g,
            "iso":        iso,
            "aperture_f": aperture,
            "exposure_s": shutter,
        }

    except Exception as e:
        print(f"FAILED — {e}")
        return None


def main():
    if not RAW_DIR.exists():
        raise FileNotFoundError(f"RAW_DIR not found:\n  {RAW_DIR}\nCheck the path at the top of this script.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cr2_files = sorted(RAW_DIR.glob("*.CR2")) + sorted(RAW_DIR.glob("*.cr2"))
    if not cr2_files:
        print(f"No CR2 files found in:\n  {RAW_DIR}")
        return

    print(f"Found {len(cr2_files)} CR2 files\n")
    results = [r for f in cr2_files if (r := process_file(f)) is not None]

    if not results:
        print("No files processed successfully.")
        return

    df_img   = pd.DataFrame(results)
    img_csv  = OUT_DIR / "sky_brightness_per_image.csv"
    df_img.to_csv(img_csv, index=False)

    # Average p1/p2/p3 per site × condition, then pivot wide
    df_by_cond = (
        df_img
        .groupby(["site", "condition"], as_index=False)
        .agg(
            mean_dsu=("mean_dsu",    "mean"),
            median_dsu=("median_dsu", "mean"),   # avg of per-image medians across p1/p2/p3
            n_images=("filename",    "count"),
        )
    )

    # Mean DSU pivot
    df_mean = df_by_cond.pivot(index="site", columns="condition", values="mean_dsu")
    df_mean.columns.name = None
    df_mean = df_mean.reset_index().rename(columns={
        "dark":     "brightness_dark",
        "white100": "brightness_white100",
    })

    # Median DSU pivot
    df_med = df_by_cond.pivot(index="site", columns="condition", values="median_dsu")
    df_med.columns.name = None
    df_med = df_med.reset_index().rename(columns={
        "dark":     "brightness_dark_median",
        "white100": "brightness_white100_median",
    })

    df_sites = df_mean.merge(df_med, on="site", how="left")

    site_csv = OUT_DIR / "sky_brightness_measurements.csv"
    df_sites.to_csv(site_csv, index=False)

    print(f"\nDone! {len(results)} images processed | {df_sites.shape[0]} sites")
    print(f"  Per-image : {img_csv}")
    print(f"  Per-site  : {site_csv}\n")
    print(df_sites.to_string(index=False))


if __name__ == "__main__":
    main()
