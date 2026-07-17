#!/usr/bin/env python3
"""Download annotated DECaLS JPEG cutouts or DS9-friendly FITS cubes."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

PIXEL_SCALE_ARCSEC = 0.262
BASE_IMAGE_SIZE = 2048
JPEG_EXPORT_SIZE_PX = 1200
SURVEY_LAYER = "ls-dr9"
DEFAULT_BRIGHTEST_COUNT = 5
DEFAULT_AUTO_DESI_RADIUS_DEG = 0.02
DEFAULT_AUTO_LS_RADIUS_DEG = 0.0166667
DOWNLOAD_RETRY_ATTEMPTS = 3
AUTO_INPUT_TOKEN = "__auto__"
DESI_OBJECT_CATALOG_FILENAME = "object_catalog.csv"
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


@dataclass
class HelperRunResult:
    """Capture the outcome of a helper-script execution."""

    success: bool
    stdout: str
    stderr: str
    returncode: int


def positive_float(value: str) -> float:
    """Parse a positive float for argparse."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than 0")
    return parsed


def ra_degrees(value: str) -> float:
    """Parse and normalize a right ascension value in degrees."""
    parsed = float(value)
    if not math.isfinite(parsed):
        raise argparse.ArgumentTypeError("right ascension must be finite")
    return parsed % 360.0


def dec_degrees(value: str) -> float:
    """Parse a declination value in degrees."""
    parsed = float(value)
    if not math.isfinite(parsed) or not -90 <= parsed <= 90:
        raise argparse.ArgumentTypeError("declination must be between -90 and 90 degrees")
    return parsed


def nonnegative_float(value: str) -> float:
    """Parse a finite, non-negative float for redshift."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed


def sanitize_name(value: str) -> str:
    """Make a user-provided label safe for filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._-") or "cluster"


def is_auto_input(path: Path | None) -> bool:
    """Return True when an overlay path flag was provided without a value."""
    return path == Path(AUTO_INPUT_TOKEN)


def build_cutout_url(ra: float, dec: float, size: int, *, fits_format: bool, pxscale: float) -> str:
    """Build the Legacy Survey cutout URL."""
    if fits_format:
        return (
            "https://www.legacysurvey.org/viewer/fits-cutout"
            f"?ra={ra}&dec={dec}&pixscale={pxscale}&layer={SURVEY_LAYER}&size={size}"
        )

    return (
        "https://www.legacysurvey.org/viewer/jpeg-cutout"
        f"?ra={ra}&dec={dec}&size={size}&layer={SURVEY_LAYER}&pixscale={pxscale}&bands=grz"
    )


def build_default_ls_catalog_path(output_dir: Path, name: str, ra: float, dec: float, radius_deg: float) -> Path:
    """Return the default LS DR10 output CSV path."""
    safe_name = sanitize_name(name)
    filename = f"ls_dr10_{safe_name}_ra{ra:.5f}_dec{dec:.5f}_r{radius_deg:.5f}.csv"
    return output_dir / filename


def run_helper_script(
    script_name: str,
    arguments: list[str],
    *,
    allow_failure: bool = False,
) -> HelperRunResult:
    """Run another repository script using the current Python interpreter."""
    script_path = Path(__file__).with_name(script_name)
    command = [sys.executable, str(script_path), *arguments]
    print(f"Running helper: {' '.join(command)}")
    try:
        completed = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        if stdout.strip():
            print(stdout.rstrip())
        if stderr.strip():
            print(stderr.rstrip(), file=sys.stderr)
        if allow_failure:
            print(f"Helper script failed, continuing without it: {script_name}")
            return HelperRunResult(
                success=False,
                stdout=stdout,
                stderr=stderr,
                returncode=exc.returncode,
            )
        raise SystemExit(f"Helper script failed: {script_name}") from exc
    if completed.stdout.strip():
        print(completed.stdout.rstrip())
    if completed.stderr.strip():
        print(completed.stderr.rstrip(), file=sys.stderr)
    return HelperRunResult(
        success=True,
        stdout=completed.stdout,
        stderr=completed.stderr,
        returncode=completed.returncode,
    )


def describe_desi_auto_result(result: HelperRunResult | None, catalog_path: Path) -> str:
    """Explain why the auto-generated DESI catalog is unavailable."""
    if result is None:
        return (
            "DESI auto catalog is unavailable because no object_catalog.csv was found at "
            f"{catalog_path}."
        )

    helper_output = "\n".join(part for part in [result.stdout, result.stderr] if part).lower()
    if not result.success and any(token in helper_output for token in ["readtimeout", "timeouterror", "timed out"]):
        return (
            "DESI auto catalog could not be created because the SPARCL/NOIRLab request timed out. "
            "This sometimes happens with the remote service. You can try the same command again."
        )

    if "no objects found." in helper_output:
        return (
            "DESI auto catalog was not created because the DESI search completed but found no matching "
            f"objects within the default {DEFAULT_AUTO_DESI_RADIUS_DEG:.2f} deg radius."
        )

    if not result.success:
        return (
            "DESI auto catalog could not be created because the DESI helper script failed. "
            "You can try the command again, or run desi_download_spectra.py manually for more detail."
        )

    return (
        "DESI auto catalog was not created even though the DESI helper finished. "
        "You can try the command again."
    )


def ensure_overlay_inputs(
    *,
    ra: float,
    dec: float,
    name: str,
    desi_csv: Path | None,
    ls10_csv: Path | None,
) -> tuple[Path | None, Path | None]:
    """Resolve overlay inputs, auto-generating missing catalogs when requested."""
    resolved_desi_csv = desi_csv
    resolved_ls10_csv = ls10_csv
    default_output_dir = Path(name).expanduser().resolve()

    if is_auto_input(desi_csv):
        default_output_dir.mkdir(parents=True, exist_ok=True)
        resolved_desi_csv = default_output_dir / DESI_OBJECT_CATALOG_FILENAME
        desi_helper_result: HelperRunResult | None = None
        if not resolved_desi_csv.exists():
            desi_helper_result = run_helper_script(
                "desi_download_spectra.py",
                [
                    "--ra",
                    str(ra),
                    "--dec",
                    str(dec),
                    "--radius",
                    str(DEFAULT_AUTO_DESI_RADIUS_DEG),
                    "--output",
                    str(default_output_dir),
                    "--overwrite",
                ],
                allow_failure=True,
            )
        if not resolved_desi_csv.exists():
            print(describe_desi_auto_result(desi_helper_result, resolved_desi_csv))
            print("Continuing without DESI markers.")
            resolved_desi_csv = None

    if is_auto_input(ls10_csv):
        default_output_dir.mkdir(parents=True, exist_ok=True)
        resolved_ls10_csv = build_default_ls_catalog_path(
            default_output_dir,
            name,
            ra,
            dec,
            DEFAULT_AUTO_LS_RADIUS_DEG,
        )
        if not resolved_ls10_csv.exists():
            run_helper_script(
                "ls_dr10_catalog_download.py",
                [
                    "--ra",
                    str(ra),
                    "--dec",
                    str(dec),
                    "--name",
                    name,
                    "--output-dir",
                    str(default_output_dir),
                ],
            )
        if not resolved_ls10_csv.exists():
            raise SystemExit(f"Expected LS catalog was not created: {resolved_ls10_csv}")

    return resolved_desi_csv, resolved_ls10_csv


def download_file(url: str, destination: Path) -> None:
    """Download a remote file atomically, with a bounded network timeout."""
    print(f"Downloading: {url}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary_path = Path(handle.name)
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Astro-Scripts/1.0"})
        for attempt in range(1, DOWNLOAD_RETRY_ATTEMPTS + 1):
            try:
                with urllib.request.urlopen(request, timeout=60) as response, temporary_path.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
                temporary_path.replace(destination)
                return
            except (urllib.error.URLError, OSError, TimeoutError):
                if attempt == DOWNLOAD_RETRY_ATTEMPTS:
                    raise
                print(f"Download failed on attempt {attempt}/{DOWNLOAD_RETRY_ATTEMPTS}; retrying...")
                time.sleep(attempt)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def compute_physical_scale_kpc(redshift: float) -> float | None:
    """Return the physical size subtended by 1 arcmin, in kpc."""
    if math.isclose(redshift, 0.0):
        return None

    scale = cosmo.angular_diameter_distance(redshift) * (1 * u.arcmin).to(u.radian).value
    return scale.to(u.kpc).value


def draw_scale_bar(*, image_size: int, redshift: float, name: str) -> None:
    """Draw a 1 arcmin scale bar and annotation block."""
    bar_arcsec = 60.0
    bar_length_px = bar_arcsec / PIXEL_SCALE_ARCSEC

    margin_x = image_size * 0.08
    margin_y = image_size * 0.10
    x_start = margin_x
    x_end = x_start + bar_length_px
    y_bar = image_size - margin_y

    plt.plot([x_start, x_end], [y_bar, y_bar], color="white", lw=2.2, solid_capstyle="butt")
    plt.plot([x_start, x_start], [y_bar - 8, y_bar + 8], color="white", lw=2.0)
    plt.plot([x_end, x_end], [y_bar - 8, y_bar + 8], color="white", lw=2.0)

    physical_scale = compute_physical_scale_kpc(redshift)
    if physical_scale is None:
        scale_text = "1 arcmin"
        redshift_text = f"{name} | z unknown"
    else:
        scale_text = f"1 arcmin = {physical_scale:.0f} kpc"
        redshift_text = f"{name} | z = {redshift:.3f}"

    plt.text(
        x_start,
        y_bar - 22,
        scale_text,
        color="white",
        fontsize=12,
        ha="left",
        va="bottom",
        bbox={"facecolor": "black", "alpha": 0.35, "pad": 2.0, "edgecolor": "none"},
    )
    plt.text(
        x_start,
        y_bar + 18,
        redshift_text,
        color="white",
        fontsize=12,
        ha="left",
        va="top",
        bbox={"facecolor": "black", "alpha": 0.35, "pad": 2.0, "edgecolor": "none"},
    )


def format_axes(fig: plt.Figure) -> None:
    """Remove axis decorations from all axes in the figure."""
    for axis in fig.axes:
        axis.tick_params(labelbottom=False, labelleft=False)
        axis.axis("off")


def make_ds9_friendly_fits(input_path: Path, output_path: Path) -> Path:
    """Reorder a DECaLS RGB FITS cube from g,r,z to z,r,g for DS9."""
    with fits.open(input_path) as hdul:
        primary_hdu = hdul[0]
        data = np.asarray(primary_hdu.data)

        if data.ndim != 3:
            raise ValueError(f"Expected a 3D FITS cube, found shape {data.shape}")
        if data.shape[0] != 3:
            raise ValueError(f"Expected 3 bands on axis 0, found shape {data.shape}")

        reordered_data = data[[2, 1, 0], :, :]
        output_hdu = fits.PrimaryHDU(data=reordered_data, header=primary_hdu.header.copy())
        output_hdu.writeto(output_path, overwrite=True)

    return output_path


def read_object_catalog(csv_path: Path) -> list[dict[str, object]]:
    """Read DESI object positions from an object_catalog.csv file."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Catalog file not found: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Catalog file is empty: {csv_path}")

        required = {"sparcl_id", "ra", "dec"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Catalog file {csv_path} is missing required columns: {sorted(missing)}"
            )

        rows: list[dict[str, object]] = []
        for row in reader:
            redshift_value = row.get("redshift")
            redshift: float | None = None
            if redshift_value not in (None, ""):
                try:
                    parsed_redshift = float(redshift_value)
                except (TypeError, ValueError):
                    parsed_redshift = math.nan
                if math.isfinite(parsed_redshift):
                    redshift = parsed_redshift

            rows.append(
                {
                    "sparcl_id": str(row["sparcl_id"]),
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "redshift": redshift,
                }
            )
    return rows


def read_ls10_objects(csv_path: Path, brightest_count: int) -> list[dict[str, object]]:
    """Read an LS DR10 CSV, sort by mag_i, and keep the brightest rows."""
    if not csv_path.exists():
        raise FileNotFoundError(f"LS catalog file not found: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"LS catalog file is empty: {csv_path}")

        required = {"ra", "dec", "mag_i"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"LS catalog file {csv_path} is missing required columns: {sorted(missing)}"
            )

        rows: list[dict[str, object]] = []
        for row in reader:
            try:
                mag_i = float(row["mag_i"])
                ra = float(row["ra"])
                dec = float(row["dec"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(mag_i):
                continue
            rows.append({"ra": ra, "dec": dec, "mag_i": mag_i})

    rows.sort(key=lambda row: row["mag_i"])
    return rows[:brightest_count]


def plot_catalog_markers(
    *,
    catalog_rows: list[dict[str, object]],
    ra_center: float,
    dec_center: float,
    image_size: int,
    label_chars: int,
) -> None:
    """Overlay object markers from a DESI catalog on the JPEG image."""
    x_center = image_size / 2
    y_center = image_size / 2
    cos_dec = math.cos(math.radians(dec_center))

    for row in catalog_rows:
        delta_ra_arcsec = (float(row["ra"]) - ra_center) * cos_dec * 3600.0
        delta_dec_arcsec = (float(row["dec"]) - dec_center) * 3600.0

        # Assume the cutout is north-up and east-left.
        x_pos = x_center - (delta_ra_arcsec / PIXEL_SCALE_ARCSEC)
        y_pos = y_center - (delta_dec_arcsec / PIXEL_SCALE_ARCSEC)

        if not (0 <= x_pos < image_size and 0 <= y_pos < image_size):
            continue

        label = str(row["sparcl_id"])[:label_chars]
        redshift = row.get("redshift")
        if isinstance(redshift, (int, float)) and math.isfinite(float(redshift)):
            label = f"{label} z={float(redshift):.2f}"
        plt.scatter(
            x_pos,
            y_pos,
            s=100,
            facecolors="none",
            edgecolors="cyan",
            linewidths=1.2,
        )
        plt.text(
            x_pos + 10,
            y_pos - 10,
            label,
            color="cyan",
            fontsize=10,
            weight="bold",
            bbox={"facecolor": "black", "alpha": 0.35, "pad": 1.5, "edgecolor": "none"},
        )


def plot_ls10_markers(
    *,
    catalog_rows: list[dict[str, object]],
    ra_center: float,
    dec_center: float,
    image_size: int,
) -> None:
    """Overlay the brightest LS DR10 objects and label them by mag_i."""
    x_center = image_size / 2
    y_center = image_size / 2
    cos_dec = math.cos(math.radians(dec_center))

    for row in catalog_rows:
        delta_ra_arcsec = (float(row["ra"]) - ra_center) * cos_dec * 3600.0
        delta_dec_arcsec = (float(row["dec"]) - dec_center) * 3600.0
        x_pos = x_center - (delta_ra_arcsec / PIXEL_SCALE_ARCSEC)
        y_pos = y_center - (delta_dec_arcsec / PIXEL_SCALE_ARCSEC)

        if not (0 <= x_pos < image_size and 0 <= y_pos < image_size):
            continue

        plt.scatter(
            x_pos,
            y_pos,
            s=110,
            facecolors="none",
            edgecolors="yellow",
            linewidths=1.6,
        )
        plt.text(
            x_pos + 10,
            y_pos + 10,
            f"{float(row['mag_i']):.1f}",
            color="yellow",
            fontsize=10,
            weight="bold",
            bbox={"facecolor": "black", "alpha": 0.4, "pad": 1.5, "edgecolor": "none"},
        )


def create_annotated_jpeg(
    *,
    temp_image_path: Path,
    output_path: Path,
    name: str,
    redshift: float,
    fraction_size: float,
    ra_center: float,
    dec_center: float,
    desi_csv: Path | None,
    ls10_csv: Path | None,
    brightest_count: int,
    label_chars: int,
) -> Path:
    """Annotate a downloaded JPEG cutout and save it."""
    image = mpimg.imread(temp_image_path)
    image_size = image.shape[0]

    figure = plt.figure(
        figsize=(JPEG_EXPORT_SIZE_PX / 120, JPEG_EXPORT_SIZE_PX / 120),
        dpi=120,
    )
    axis = figure.add_axes([0, 0, 1, 1])
    axis.imshow(image)

    if desi_csv is not None:
        catalog_rows = read_object_catalog(desi_csv)
        plot_catalog_markers(
            catalog_rows=catalog_rows,
            ra_center=ra_center,
            dec_center=dec_center,
            image_size=image_size,
            label_chars=label_chars,
        )

    if ls10_csv is not None:
        ls10_rows = read_ls10_objects(ls10_csv, brightest_count)
        plot_ls10_markers(
            catalog_rows=ls10_rows,
            ra_center=ra_center,
            dec_center=dec_center,
            image_size=image_size,
        )

    draw_scale_bar(image_size=image_size, redshift=redshift, name=name)

    format_axes(figure)
    plt.savefig(output_path, dpi=120)
    plt.close()
    return output_path


def build_jpeg_output_path(
    *,
    index: int,
    name: str,
    desi_csv: Path | None,
    ls10_csv: Path | None,
) -> Path:
    """Build a JPEG output path that keeps overlay variants separate."""
    suffix_parts: list[str] = []
    if desi_csv is not None:
        suffix_parts.append("desi")
    if ls10_csv is not None:
        suffix_parts.append("lsbright")

    suffix = ""
    if suffix_parts:
        suffix = "_" + "_".join(suffix_parts)

    return Path(f"img_ix{index:05d}_annoted_{name}{suffix}.jpg")


def download_decal_image(
    index: int,
    ra: float,
    dec: float,
    name: str,
    redshift: float,
    fraction_size: float = 1,
    *,
    fits_format: bool = False,
    keep_raw_fits: bool = False,
    desi_csv: Path | None = None,
    ls10_csv: Path | None = None,
    brightest_count: int = DEFAULT_BRIGHTEST_COUNT,
    label_chars: int = 3,
    overwrite: bool = False,
) -> Path:
    """Download a DECaLS cutout as a JPEG or a DS9-friendly FITS file."""
    size = int(BASE_IMAGE_SIZE / fraction_size)
    if size < 1:
        raise ValueError(f"fraction_size must be no greater than {BASE_IMAGE_SIZE}")
    url = build_cutout_url(ra, dec, size, fits_format=fits_format, pxscale=PIXEL_SCALE_ARCSEC)

    if fits_format:
        output_path = Path(f"img_ix{index:05d}_annoted_{name}.fits")
        raw_path = output_path.with_name(f"{output_path.stem}_raw.fits")
        if output_path.exists() and not overwrite:
            print(f"Output already exists, skipping download: {output_path}")
            return output_path

        try:
            download_file(url, raw_path)
            make_ds9_friendly_fits(raw_path, output_path)
        except (urllib.error.URLError, OSError, ValueError) as exc:
            raise SystemExit(f"Failed to create FITS cutout: {exc}") from exc
        finally:
            if raw_path.exists() and not keep_raw_fits:
                raw_path.unlink()

        print(f"Output filename is {output_path}")
        if keep_raw_fits:
            print(f"Kept original downloaded FITS as {raw_path}")
        return output_path

    output_path = build_jpeg_output_path(
        index=index,
        name=name,
        desi_csv=desi_csv,
        ls10_csv=ls10_csv,
    )
    if output_path.exists() and not overwrite:
        print(f"Output already exists, skipping download: {output_path}")
        return output_path

    temp_file = tempfile.NamedTemporaryFile(suffix=".jpg", delete=False)
    temp_path = Path(temp_file.name)
    temp_file.close()

    try:
        download_file(url, temp_path)
        create_annotated_jpeg(
            temp_image_path=temp_path,
            output_path=output_path,
            name=name,
            redshift=redshift,
            fraction_size=fraction_size,
            ra_center=ra,
            dec_center=dec,
            desi_csv=desi_csv,
            ls10_csv=ls10_csv,
            brightest_count=brightest_count,
            label_chars=label_chars,
        )
    except (urllib.error.URLError, OSError, ValueError, FileNotFoundError) as exc:
        raise SystemExit(f"Failed to create JPEG cutout: {exc}") from exc
    finally:
        if temp_path.exists():
            temp_path.unlink()

    print(f"Output filename is {output_path}")
    return output_path


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Download DECaLS cutouts as annotated JPEG files or FITS cubes with "
            "band ordering fixed for DS9."
        )
    )
    parser.add_argument("index", type=int, help="Identifier index or row number.")
    parser.add_argument("ra", type=ra_degrees, help="Right ascension in degrees.")
    parser.add_argument("dec", type=dec_degrees, help="Declination in degrees.")
    parser.add_argument("name", help="Object name used in the output filename.")
    parser.add_argument(
        "redshift",
        type=nonnegative_float,
        nargs="?",
        default=0.2,
        help="Target redshift used for the physical scale label (default: 0.2). Use 0 if it is unknown.",
    )
    parser.add_argument(
        "fraction_size",
        type=positive_float,
        nargs="?",
        default=1,
        help="Image scale factor: 1=2048 px, 2=1024 px, 4=512 px (default: 1).",
    )
    parser.add_argument(
        "--fits",
        action="store_true",
        help="Download a FITS cube instead of an annotated JPEG.",
    )
    parser.add_argument(
        "--keep-raw-fits",
        action="store_true",
        help="Keep the original downloaded FITS cube alongside the DS9-friendly version.",
    )
    parser.add_argument(
        "--desi-csv",
        dest="desi_csv",
        nargs="?",
        type=Path,
        const=Path(AUTO_INPUT_TOKEN),
        help=(
            "Path to object_catalog.csv from desi_download_spectra.py. "
            "If you provide the flag without a value, the script auto-creates the default DESI catalog in <name>/."
        ),
    )
    parser.add_argument(
        "--ls10-csv",
        dest="ls10_csv",
        nargs="?",
        type=Path,
        const=Path(AUTO_INPUT_TOKEN),
        help=(
            "CSV file produced by ls_dr10_catalog_download.py. "
            "If you provide the flag without a value, the script auto-creates the default LS catalog in <name>/."
        ),
    )
    parser.add_argument(
        "--brightest-count",
        type=int,
        default=DEFAULT_BRIGHTEST_COUNT,
        help=f"Number of brightest LS objects to mark from --ls10-csv (default: {DEFAULT_BRIGHTEST_COUNT}).",
    )
    parser.add_argument(
        "--label-chars",
        type=int,
        default=3,
        help="Number of leading sparcl_id characters to show beside each marker (default: 3).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recreate the final JPEG even if it already exists.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.label_chars <= 0:
        raise SystemExit("--label-chars must be greater than 0")
    if args.brightest_count <= 0:
        raise SystemExit("--brightest-count must be greater than 0")

    desi_csv = args.desi_csv.expanduser().resolve() if args.desi_csv and not is_auto_input(args.desi_csv) else args.desi_csv
    ls10_csv = args.ls10_csv.expanduser().resolve() if args.ls10_csv and not is_auto_input(args.ls10_csv) else args.ls10_csv
    desi_csv, ls10_csv = ensure_overlay_inputs(
        ra=args.ra,
        dec=args.dec,
        name=args.name,
        desi_csv=desi_csv,
        ls10_csv=ls10_csv,
    )

    if not args.fits:
        planned_output_path = build_jpeg_output_path(
            index=args.index,
            name=args.name,
            desi_csv=desi_csv,
            ls10_csv=ls10_csv,
        )
        if planned_output_path.exists() and not args.overwrite:
            print(f"Output already exists, skipping download: {planned_output_path}")
            return

    download_decal_image(
        args.index,
        args.ra,
        args.dec,
        args.name,
        redshift=args.redshift,
        fraction_size=args.fraction_size,
        fits_format=args.fits,
        keep_raw_fits=args.keep_raw_fits,
        desi_csv=desi_csv,
        ls10_csv=ls10_csv,
        brightest_count=args.brightest_count,
        label_chars=args.label_chars,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
