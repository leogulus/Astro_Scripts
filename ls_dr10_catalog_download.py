#!/usr/bin/env python3
"""Download Legacy Surveys DR10 Tractor catalog data for a sky box."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

from astroquery.utils.tap.core import TapPlus

TAP_URL = "https://datalab.noirlab.edu/tap"
DEFAULT_RADIUS_DEG = 0.0166667
DEFAULT_OUTPUT_DIR = Path(".")


def positive_float(value: str) -> float:
    """Parse a positive float for argparse."""
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than 0")
    return parsed


def dec_degrees(value: str) -> float:
    """Parse a declination value in degrees."""
    parsed = float(value)
    if not -90 <= parsed <= 90:
        raise argparse.ArgumentTypeError("declination must be between -90 and 90 degrees")
    return parsed


def sanitize_name(value: str) -> str:
    """Make a user-provided label safe for filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._-") or "cluster"


def build_query(ra_min: float, ra_max: float, dec_min: float, dec_max: float) -> str:
    """Build a fast RA/Dec box query for LS DR10."""
    return f"""
    SELECT
        ra, dec,
        flux_g, flux_r, flux_i, flux_z,
        flux_w1, flux_w2,
        mag_g, mag_r, mag_i, mag_z,
        mag_w1, mag_w2,
        dered_flux_g, dered_flux_r, dered_flux_i, dered_flux_z,
        dered_flux_w1, dered_flux_w2,
        dered_mag_g, dered_mag_r, dered_mag_i, dered_mag_z,
        dered_mag_w1, dered_mag_w2,
        flux_ivar_g, flux_ivar_r, flux_ivar_i, flux_ivar_z,
        flux_ivar_w1, flux_ivar_w2,
        snr_g, snr_r, snr_i, snr_z,
        snr_w1, snr_w2,
        mw_transmission_g, mw_transmission_r, mw_transmission_i,
        mw_transmission_z, mw_transmission_w1, mw_transmission_w2,
        type, maskbits, shape_r,
        fracflux_g, fracflux_r, fracflux_i, fracflux_z
    FROM ls_dr10.tractor
    WHERE ra BETWEEN {ra_min:.10f} AND {ra_max:.10f}
      AND dec BETWEEN {dec_min:.10f} AND {dec_max:.10f}
    """


def compute_box_bounds(ra_deg: float, dec_deg: float, radius_deg: float) -> tuple[float, float, float, float]:
    """Return the RA/Dec bounds for a box centered on the target position."""
    dec_min = dec_deg - radius_deg
    dec_max = dec_deg + radius_deg

    cos_dec = max(abs(math.cos(math.radians(dec_deg))), 1e-6)
    ra_delta = radius_deg / cos_dec
    ra_min = ra_deg - ra_delta
    ra_max = ra_deg + ra_delta

    return ra_min, ra_max, dec_min, dec_max


def build_output_path(output_dir: Path, name: str, ra_deg: float, dec_deg: float, radius_deg: float) -> Path:
    """Create a predictable output filename for the result table."""
    safe_name = sanitize_name(name)
    filename = f"ls_dr10_{safe_name}_ra{ra_deg:.5f}_dec{dec_deg:.5f}_r{radius_deg:.5f}.csv"
    return output_dir / filename


def fetch_catalog(query: str):
    """Run the TAP query and return the result table."""
    tap = TapPlus(url=TAP_URL)
    job = tap.launch_job(query)
    return job.get_results()


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Download Legacy Surveys DR10 Tractor photometry inside a fast RA/Dec box "
            "centered on the target position."
        )
    )
    parser.add_argument(
        "--name",
        type=str,
        default="cluster",
        help="Object or field name used in the output filename.",
    )
    parser.add_argument(
        "--ra",
        type=float,
        required=True,
        help="Central right ascension in degrees.",
    )
    parser.add_argument(
        "--dec",
        type=dec_degrees,
        required=True,
        help="Central declination in degrees.",
    )
    parser.add_argument(
        "--radius",
        type=positive_float,
        default=DEFAULT_RADIUS_DEG,
        help=f"Half-width of the search box in degrees (default: {DEFAULT_RADIUS_DEG}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the CSV file will be written (default: current directory).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ra_min, ra_max, dec_min, dec_max = compute_box_bounds(args.ra, args.dec, args.radius)
    query = build_query(ra_min, ra_max, dec_min, dec_max)

    print("Submitting box query to LS DR10 TAP service...")
    print(
        "Search bounds: "
        f"RA [{ra_min:.5f}, {ra_max:.5f}] deg, "
        f"Dec [{dec_min:.5f}, {dec_max:.5f}] deg"
    )

    try:
        result_table = fetch_catalog(query)
        dataframe = result_table.to_pandas()
    except Exception as exc:
        raise SystemExit(f"Query failed: {exc}") from exc

    output_path = build_output_path(output_dir, args.name, args.ra, args.dec, args.radius)
    dataframe.to_csv(output_path, index=False)

    print(f"Saved file: {output_path}")
    print(f"Number of sources: {len(dataframe)}")


if __name__ == "__main__":
    main()
