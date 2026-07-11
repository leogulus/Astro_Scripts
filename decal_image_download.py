#!/usr/bin/env python3
"""Download annotated DECaLS JPEG cutouts or DS9-friendly FITS cubes."""

from __future__ import annotations

import argparse
import csv
import math
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

PIXEL_SCALE_ARCSEC = 0.262
BASE_IMAGE_SIZE = 2048
SURVEY_LAYER = "ls-dr9"
OBJECT_CATALOG_FILENAME = "object_catalog.csv"
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def positive_float(value: str) -> float:
    """Parse a positive float for argparse."""
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than 0")
    return parsed


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


def download_file(url: str, destination: Path) -> None:
    """Download a remote file to disk."""
    print(f"Downloading: {url}")
    urllib.request.urlretrieve(url, destination)


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


def read_object_catalog(catalog_dir: Path) -> list[dict[str, object]]:
    """Read DESI object positions from object_catalog.csv."""
    catalog_path = catalog_dir / OBJECT_CATALOG_FILENAME
    if not catalog_path.exists():
        raise FileNotFoundError(f"Catalog file not found: {catalog_path}")

    with catalog_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Catalog file is empty: {catalog_path}")

        required = {"sparcl_id", "ra", "dec"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Catalog file {catalog_path} is missing required columns: {sorted(missing)}"
            )

        rows: list[dict[str, object]] = []
        for row in reader:
            rows.append(
                {
                    "sparcl_id": str(row["sparcl_id"]),
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                }
            )
    return rows


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


def create_annotated_jpeg(
    *,
    temp_image_path: Path,
    output_path: Path,
    name: str,
    redshift: float,
    fraction_size: float,
    ra_center: float,
    dec_center: float,
    catalog_dir: Path | None,
    label_chars: int,
) -> Path:
    """Annotate a downloaded JPEG cutout and save it."""
    image = mpimg.imread(temp_image_path)
    image_size = image.shape[0]

    figure = plt.figure(figsize=(10, 10))
    plt.imshow(image)

    if catalog_dir is not None:
        catalog_rows = read_object_catalog(catalog_dir)
        plot_catalog_markers(
            catalog_rows=catalog_rows,
            ra_center=ra_center,
            dec_center=dec_center,
            image_size=image_size,
            label_chars=label_chars,
        )

    draw_scale_bar(image_size=image_size, redshift=redshift, name=name)

    format_axes(figure)
    plt.tight_layout()
    plt.savefig(output_path, dpi=120)
    plt.close()
    return output_path


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
    catalog_dir: Path | None = None,
    label_chars: int = 3,
) -> Path:
    """Download a DECaLS cutout as a JPEG or a DS9-friendly FITS file."""
    size = int(BASE_IMAGE_SIZE / fraction_size)
    url = build_cutout_url(ra, dec, size, fits_format=fits_format, pxscale=PIXEL_SCALE_ARCSEC)

    if fits_format:
        output_path = Path(f"img_ix{index:05d}_annoted_{name}.fits")
        raw_path = output_path.with_name(f"{output_path.stem}_raw.fits")

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

    output_path = Path(f"img_ix{index:05d}_annoted_{name}.jpg")
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
            catalog_dir=catalog_dir,
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
    parser.add_argument("ra", type=float, help="Right ascension in degrees.")
    parser.add_argument("dec", type=float, help="Declination in degrees.")
    parser.add_argument("name", help="Object name used in the output filename.")
    parser.add_argument(
        "redshift",
        type=float,
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
        "--catalog-dir",
        type=Path,
        help=(
            "Directory containing object_catalog.csv from desi_download_spectra.py. "
            "When provided, the JPEG output marks all catalog objects using their sky offsets from the image center."
        ),
    )
    parser.add_argument(
        "--label-chars",
        type=int,
        default=3,
        help="Number of leading sparcl_id characters to show beside each marker (default: 3).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.label_chars <= 0:
        raise SystemExit("--label-chars must be greater than 0")

    download_decal_image(
        args.index,
        args.ra,
        args.dec,
        args.name,
        redshift=args.redshift,
        fraction_size=args.fraction_size,
        fits_format=args.fits,
        keep_raw_fits=args.keep_raw_fits,
        catalog_dir=args.catalog_dir.expanduser().resolve() if args.catalog_dir else None,
        label_chars=args.label_chars,
    )


if __name__ == "__main__":
    main()
