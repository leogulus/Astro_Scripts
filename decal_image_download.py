#!/usr/bin/env python3
"""Download annotated DECaLS JPEG cutouts or DS9-friendly FITS cubes."""

from __future__ import annotations

import argparse
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


def create_annotated_jpeg(
    *,
    temp_image_path: Path,
    output_path: Path,
    name: str,
    redshift: float,
    fraction_size: float,
) -> Path:
    """Annotate a downloaded JPEG cutout and save it."""
    image = mpimg.imread(temp_image_path)
    figure = plt.figure(figsize=(10, 10))
    plt.imshow(image)

    fraction_half = fraction_size / 2
    x_pos = 50 / fraction_half
    y_pos = 140 / fraction_half

    plt.annotate("1 arcmin", (x_pos + 40, y_pos - 10 / fraction_half), color="white", fontsize=13)

    physical_scale = compute_physical_scale_kpc(redshift)
    if physical_scale is None:
        plt.annotate(name, (x_pos, y_pos - 50 / fraction_half), color="white", fontsize=13)
        scale_label = "z unknown"
    else:
        scale_label = f"{physical_scale:.0f} kpc"

    plt.annotate(scale_label, (x_pos + 40, y_pos + 17 / fraction_half), color="white", fontsize=13)
    plt.annotate(
        "",
        xy=(x_pos, y_pos + 50),
        xytext=(x_pos + 229, y_pos + 50),
        arrowprops={"arrowstyle": "-", "color": "white"},
    )

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
        )
    except (urllib.error.URLError, OSError, ValueError) as exc:
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
        help="Target redshift. Use 0 if it is unknown.",
    )
    parser.add_argument(
        "fraction_size",
        type=positive_float,
        help="Image scale factor: 1=2048 px, 2=1024 px, 4=512 px.",
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    download_decal_image(
        args.index,
        args.ra,
        args.dec,
        args.name,
        redshift=args.redshift,
        fraction_size=args.fraction_size,
        fits_format=args.fits,
        keep_raw_fits=args.keep_raw_fits,
    )


if __name__ == "__main__":
    main()
