#!/usr/bin/env python3
"""Download Legacy Surveys DR10 Tractor catalog data for a sky box."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from astroquery.utils.tap.core import TapPlus

TAP_URL = "https://datalab.noirlab.edu/tap"
TRACTOR_TABLE = "ls_dr10.tractor"
DEFAULT_RADIUS_DEG = 0.0166667
DEFAULT_OUTPUT_DIR = Path(".")
DEFAULT_COLUMNS = [
    "ra",
    "dec",
    "flux_g",
    "flux_r",
    "flux_i",
    "flux_z",
    "flux_w1",
    "flux_w2",
    "mag_g",
    "mag_r",
    "mag_i",
    "mag_z",
    "mag_w1",
    "mag_w2",
    "dered_flux_g",
    "dered_flux_r",
    "dered_flux_i",
    "dered_flux_z",
    "dered_flux_w1",
    "dered_flux_w2",
    "dered_mag_g",
    "dered_mag_r",
    "dered_mag_i",
    "dered_mag_z",
    "dered_mag_w1",
    "dered_mag_w2",
    "flux_ivar_g",
    "flux_ivar_r",
    "flux_ivar_i",
    "flux_ivar_z",
    "flux_ivar_w1",
    "flux_ivar_w2",
    "snr_g",
    "snr_r",
    "snr_i",
    "snr_z",
    "snr_w1",
    "snr_w2",
    "mw_transmission_g",
    "mw_transmission_r",
    "mw_transmission_i",
    "mw_transmission_z",
    "mw_transmission_w1",
    "mw_transmission_w2",
    "type",
    "maskbits",
    "shape_r",
    "fracflux_g",
    "fracflux_r",
    "fracflux_i",
    "fracflux_z",
]


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


def sanitize_name(value: str) -> str:
    """Make a user-provided label safe for filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._-") or "cluster"


def parse_column_list(value: str) -> list[str]:
    """Parse a comma-separated list of SQL column names."""
    columns = []
    for raw_column in value.split(","):
        column = raw_column.strip()
        if not column:
            continue
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", column):
            raise argparse.ArgumentTypeError(f"invalid column name: {raw_column!r}")
        columns.append(column)
    if not columns:
        raise argparse.ArgumentTypeError("provide at least one column name")
    return columns


def dedupe_columns(columns: list[str]) -> list[str]:
    """Keep the first occurrence of each column while preserving order."""
    seen: set[str] = set()
    deduped = []
    for column in columns:
        if column not in seen:
            seen.add(column)
            deduped.append(column)
    return deduped


def build_query(
    ra_min: float,
    ra_max: float,
    dec_min: float,
    dec_max: float,
    selected_columns: list[str],
) -> str:
    """Build a fast RA/Dec box query for LS DR10."""
    column_lines = ",\n        ".join(selected_columns)
    ra_clause = build_ra_clause(ra_min, ra_max)
    return f"""
    SELECT
        {column_lines}
    FROM {TRACTOR_TABLE}
    WHERE {ra_clause}
      AND dec BETWEEN {dec_min:.10f} AND {dec_max:.10f}
    """


def build_ra_clause(ra_min: float, ra_max: float) -> str:
    """Build an RA predicate that remains correct across the 0/360-degree seam."""
    if ra_max - ra_min >= 360:
        return "1 = 1"
    normalized_min = ra_min % 360.0
    normalized_max = ra_max % 360.0
    if normalized_min <= normalized_max and 0 <= ra_min and ra_max < 360:
        return f"ra BETWEEN {normalized_min:.10f} AND {normalized_max:.10f}"
    return f"(ra >= {normalized_min:.10f} OR ra <= {normalized_max:.10f})"


def compute_box_bounds(ra_deg: float, dec_deg: float, radius_deg: float) -> tuple[float, float, float, float]:
    """Return the RA/Dec bounds for a box centered on the target position."""
    dec_min = max(-90.0, dec_deg - radius_deg)
    dec_max = min(90.0, dec_deg + radius_deg)

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


def create_tap_client() -> "TapPlus":
    """Create a TAP client for the NOIRLab service."""
    try:
        from astroquery.utils.tap.core import TapPlus
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "astroquery is required for LS DR10 TAP queries. Install it with: pip install astroquery"
        ) from exc
    return TapPlus(url=TAP_URL)


def run_query(tap: "TapPlus", query: str):
    """Run a TAP query and return the result table."""
    job = tap.launch_job(query)
    return job.get_results()


def fetch_available_columns(tap: "TapPlus") -> list[dict[str, str]]:
    """Return column metadata for the LS DR10 Tractor table."""
    metadata_query = f"""
    SELECT column_name, datatype, unit, description
    FROM TAP_SCHEMA.columns
    WHERE table_name = '{TRACTOR_TABLE}'
    ORDER BY column_name
    """
    table = run_query(tap, metadata_query)
    columns = []
    for row in table:
        columns.append(
            {
                "column_name": str(row["column_name"]),
                "datatype": "" if row["datatype"] is None else str(row["datatype"]),
                "unit": "" if row["unit"] is None else str(row["unit"]),
                "description": "" if row["description"] is None else str(row["description"]),
            }
        )
    return columns


def print_available_columns(column_metadata: list[dict[str, str]]) -> None:
    """Print the available LS DR10 Tractor columns in a compact table."""
    if not column_metadata:
        print("No columns were returned by TAP_SCHEMA.columns.")
        return

    name_width = max(len(item["column_name"]) for item in column_metadata)
    type_width = max(len(item["datatype"]) for item in column_metadata)
    unit_width = max(len(item["unit"]) for item in column_metadata)

    header = (
        f"{'column_name'.ljust(name_width)}  "
        f"{'datatype'.ljust(type_width)}  "
        f"{'unit'.ljust(unit_width)}  description"
    )
    print(header)
    print("-" * len(header))

    for item in column_metadata:
        print(
            f"{item['column_name'].ljust(name_width)}  "
            f"{item['datatype'].ljust(type_width)}  "
            f"{item['unit'].ljust(unit_width)}  "
            f"{item['description']}"
        )


def validate_requested_columns(
    requested_columns: list[str],
    column_metadata: list[dict[str, str]],
) -> None:
    """Ensure every requested column exists in the remote Tractor table."""
    available = {item["column_name"] for item in column_metadata}
    missing = [column for column in requested_columns if column not in available]
    if missing:
        raise SystemExit(
            "Unknown column(s) requested: "
            f"{', '.join(missing)}. "
            "Use --list-columns to inspect the currently available ls_dr10.tractor columns."
        )


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
        type=ra_degrees,
        help="Central right ascension in degrees.",
    )
    parser.add_argument(
        "--dec",
        type=dec_degrees,
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
    parser.add_argument(
        "--extra-columns",
        type=parse_column_list,
        default=[],
        help=(
            "Comma-separated list of additional ls_dr10.tractor columns to append to the default export. "
            "Use --list-columns to inspect available options."
        ),
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Query TAP_SCHEMA.columns and print the currently available ls_dr10.tractor columns, then exit.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    tap = create_tap_client()

    if args.list_columns:
        try:
            print_available_columns(fetch_available_columns(tap))
        except Exception as exc:
            raise SystemExit(f"Could not fetch available columns: {exc}") from exc
        return

    if args.ra is None or args.dec is None:
        raise SystemExit("the following arguments are required for catalog downloads: --ra, --dec")

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_columns = dedupe_columns(DEFAULT_COLUMNS + args.extra_columns)

    if args.extra_columns:
        try:
            available_columns = fetch_available_columns(tap)
            validate_requested_columns(selected_columns, available_columns)
        except Exception as exc:
            if isinstance(exc, SystemExit):
                raise
            raise SystemExit(f"Could not validate requested columns: {exc}") from exc

    ra_min, ra_max, dec_min, dec_max = compute_box_bounds(args.ra, args.dec, args.radius)
    query = build_query(ra_min, ra_max, dec_min, dec_max, selected_columns)

    print("Submitting box query to LS DR10 TAP service...")
    print(
        "Search bounds: "
        f"RA [{ra_min:.5f}, {ra_max:.5f}] deg, "
        f"Dec [{dec_min:.5f}, {dec_max:.5f}] deg"
    )
    if args.extra_columns:
        print(f"Including extra columns: {', '.join(args.extra_columns)}")

    try:
        result_table = run_query(tap, query)
        dataframe = result_table.to_pandas()
    except Exception as exc:
        raise SystemExit(f"Query failed: {exc}") from exc

    output_path = build_output_path(output_dir, args.name, args.ra, args.dec, args.radius)
    dataframe.to_csv(output_path, index=False)

    print(f"Saved file: {output_path}")
    print(f"Number of sources: {len(dataframe)}")


if __name__ == "__main__":
    main()
