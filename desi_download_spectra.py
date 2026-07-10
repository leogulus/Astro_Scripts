#!/usr/bin/env python3
"""Download DESI/SPARCL spectra for one target or a batch of targets."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve

DEFAULT_DATA_RELEASE = "DESI-DR1"
DEFAULT_LIMIT = 1000
DEFAULT_OUTPUT_DIR = "desi_output"
DEFAULT_SMOOTH_SIGMA = 5.0
SENTINEL_FILENAME = "download_complete.txt"
CATALOG_FILENAME = "object_catalog.csv"


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


def metadata_value(meta, key: str, index: int):
    """Return a metadata value for one spectrum, handling array-like entries."""
    value = meta.get(key)
    if value is None:
        return None
    try:
        return value[index]
    except (TypeError, KeyError, IndexError):
        return value


def default_lines(include_absorption: bool) -> dict[str, float]:
    """Return rest-frame spectral lines to overlay on plots."""
    lines = {
        "[O II]": 3727.0,
        "Hbeta": 4861.33,
        "[O III]": 5006.84,
        "Halpha": 6562.80,
        "[N II]": 6583.45,
        "[S II]": 6716.44,
    }

    if include_absorption:
        lines.update(
            {
                "Ca K": 3933.66,
                "Ca H": 3968.47,
                "Hdelta": 4101.74,
                "Hgamma": 4340.47,
                "Mg b": 5175.27,
                "Na D": 5895.92,
            }
        )

    return lines


def cone_query_sql(
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    table: str,
    extra_where: str,
    limit: int,
) -> str:
    """Build the SPARCL cone-search SQL query."""
    cols = [
        "sparcl_id",
        "specid",
        "ra",
        "dec",
        "redshift",
        "redshift_err",
        "spectype",
        "data_release",
    ]

    cos_dec = max(np.cos(np.radians(dec_deg)), 1e-6)
    ra_width = radius_deg / cos_dec
    ra_min = ra_deg - ra_width
    ra_max = ra_deg + ra_width
    dec_min = dec_deg - radius_deg
    dec_max = dec_deg + radius_deg

    return f"""
    SELECT {", ".join(cols)}
    FROM {table}
    WHERE
        ra BETWEEN {ra_min:.10f} AND {ra_max:.10f}
        AND dec BETWEEN {dec_min:.10f} AND {dec_max:.10f}
        AND acos(
            sin(radians(dec))*sin(radians({dec_deg:.10f})) +
            cos(radians(dec))*cos(radians({dec_deg:.10f}))*
            cos(radians(ra-{ra_deg:.10f}))
        ) <= radians({radius_deg:.10f})
        AND {extra_where}
    LIMIT {limit}
    """


def cone_query(
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    *,
    table: str,
    data_release: str,
    limit: int,
):
    """Run the SPARCL cone search and return a pandas DataFrame."""
    from dl import queryClient as qc

    sql = cone_query_sql(
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        radius_deg=radius_deg,
        table=table,
        extra_where=f"data_release='{data_release}'",
        limit=limit,
    )
    return qc.query(sql=sql, fmt="pandas")


def create_sparcl_client():
    """Create the SPARCL API client lazily."""
    from sparcl.client import SparclClient

    return SparclClient()


def safe_float(value, default: float = 0.0) -> float:
    """Convert a metadata value to float when possible."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def build_catalog_lookup(found) -> dict[str, dict[str, object]]:
    """Build a lookup of catalog rows keyed by SPARCL ID."""
    lookup: dict[str, dict[str, object]] = {}
    for row in found.to_dict(orient="records"):
        sparcl_id = str(row["sparcl_id"])
        lookup[sparcl_id] = row
    return lookup


def save_spectrum(results, index: int, output_dir: Path) -> Path:
    """Save one spectrum as a compressed NumPy archive."""
    sparcl_id = results.meta["sparcl_id"][index]

    def array_meta(key: str):
        value = results.meta.get(key)
        if value is None:
            return None
        return np.asarray(value[index])

    output_path = output_dir / f"spectrum_{sparcl_id}.npz"
    np.savez(
        output_path,
        wavelength=np.asarray(results.spectral_axis),
        flux=np.asarray(results.flux[index]),
        model=array_meta("model"),
        ivar=array_meta("ivar"),
        mask=array_meta("mask"),
        wave_sigma=array_meta("wave_sigma"),
        sparcl_id=sparcl_id,
        redshift=metadata_value(results.meta, "redshift", index),
        specid=metadata_value(results.meta, "specid", index),
        spectype=metadata_value(results.meta, "spectype", index),
        ra=metadata_value(results.meta, "ra", index),
        dec=metadata_value(results.meta, "dec", index),
    )
    return output_path


def plot_spectrum(
    results,
    index: int,
    output_dir: Path,
    catalog_lookup: dict[str, dict[str, object]],
    *,
    show_model: bool,
    show_smooth: bool,
    smooth_sigma: float,
    include_absorption: bool,
) -> Path:
    """Plot one downloaded spectrum and save it to disk."""
    sparcl_id = results.meta["sparcl_id"][index]
    output_path = output_dir / f"spectrum_{sparcl_id}.png"

    wavelength = np.asarray(results.spectral_axis)
    flux = np.asarray(results.flux[index]) if np.ndim(results.flux) != 1 else np.asarray(results.flux)
    wave_unit = getattr(results.spectral_axis, "unit", "")
    flux_unit = getattr(results.flux, "unit", "")

    metadata = {key: metadata_value(results.meta, key, index) for key in results.meta}
    catalog_row = catalog_lookup.get(str(sparcl_id), {})
    redshift = safe_float(metadata.get("redshift"), default=0.0)
    ra = safe_float(metadata.get("ra"), default=safe_float(catalog_row.get("ra"), default=float("nan")))
    dec = safe_float(metadata.get("dec"), default=safe_float(catalog_row.get("dec"), default=float("nan")))

    plt.figure(figsize=(10, 5))
    plt.plot(wavelength, flux, alpha=0.3, lw=0.8, label="Observed")

    if show_smooth:
        kernel = Gaussian1DKernel(smooth_sigma)
        smoothed = convolve(flux, kernel, boundary="extend")
        plt.plot(wavelength, smoothed, lw=1.2, label="Smoothed")

    if show_model and "model" in results.meta:
        model = np.asarray(results.meta["model"][index])
        if model.ndim > 1:
            model = model[0]
        plt.plot(wavelength, model, lw=1.5, label="Model")

    ymax = plt.ylim()[1]
    for line_name, rest_wavelength in default_lines(include_absorption).items():
        observed_wavelength = rest_wavelength * (1 + redshift)
        plt.axvline(observed_wavelength, linestyle="--", alpha=0.3)
        plt.text(
            observed_wavelength,
            ymax * 0.95,
            line_name,
            rotation=90,
            fontsize=8,
            ha="right",
            va="top",
        )

    wave_label = f"Wavelength [{wave_unit}]" if wave_unit else "Wavelength"
    flux_label = f"Flux [{flux_unit}]" if flux_unit else "Flux"
    plt.xlabel(wave_label)
    plt.ylabel(flux_label)

    title_ra = f"{ra:.5f}" if math.isfinite(ra) else "unknown"
    title_dec = f"{dec:.5f}" if math.isfinite(dec) else "unknown"
    plt.title(
        f"SPARCL ID = {metadata.get('sparcl_id', sparcl_id)}\n"
        f"z = {redshift:.4f}, RA = {title_ra}, Dec = {title_dec}"
    )

    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()
    return output_path


def read_targets_csv(csv_path: Path) -> list[dict[str, str]]:
    """Read batch targets from a CSV file."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"CSV file not found: {csv_path}\n"
            "Provide a valid path to a CSV with columns: name, ra, dec, radius."
        )

    if not csv_path.is_file():
        raise ValueError(f"CSV path is not a file: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV file is empty or missing a header row")

        required = {"name", "ra", "dec", "radius"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV file is missing required columns: {sorted(missing)}")

        return list(reader)


def include_fields(show_model: bool) -> list[str]:
    """Return SPARCL fields to request for retrieval."""
    fields = [
        "sparcl_id",
        "specid",
        "data_release",
        "redshift",
        "flux",
        "wavelength",
        "ivar",
        "mask",
        "spectype",
        "ra",
        "dec",
        "wave_sigma",
    ]
    if show_model:
        fields.append("model")
    return fields


def process_target(
    *,
    client,
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    output_dir: Path,
    target_label: str,
    target_index: int,
    total_targets: int,
    table: str,
    data_release: str,
    limit: int,
    overwrite: bool,
    save_plots: bool,
    show_model: bool,
    show_smooth: bool,
    smooth_sigma: float,
    include_absorption: bool,
) -> None:
    """Download all spectra for one target and save outputs."""
    sentinel_path = output_dir / SENTINEL_FILENAME

    if output_dir.exists() and sentinel_path.exists() and not overwrite:
        print(f"Skipping {output_dir}: already completed.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Target {target_index}/{total_targets}: {target_label} ===")
    print(f"RA={ra_deg:.6f}, Dec={dec_deg:.6f}, Radius={radius_deg:.4f}")

    try:
        found = cone_query(
            ra_deg,
            dec_deg,
            radius_deg,
            table=table,
            data_release=data_release,
            limit=limit,
        )
    except Exception as exc:
        raise SystemExit(f"Cone search failed for {target_label}: {exc}") from exc

    if len(found) == 0:
        print("  No objects found.")
        return

    catalog_path = output_dir / CATALOG_FILENAME
    found.to_csv(catalog_path, index=False)

    sparcl_ids = list(found["sparcl_id"])
    catalog_lookup = build_catalog_lookup(found)
    try:
        results = client.retrieve(
            uuid_list=sparcl_ids,
            include=include_fields(show_model),
            fmt="specutils",
        )
    except Exception as exc:
        raise SystemExit(f"Spectrum retrieval failed for {target_label}: {exc}") from exc

    for index, sparcl_id in enumerate(sparcl_ids, start=1):
        save_spectrum(results, index - 1, output_dir)
        if save_plots:
            plot_spectrum(
                results,
                index - 1,
                output_dir,
                catalog_lookup,
                show_model=show_model,
                show_smooth=show_smooth,
                smooth_sigma=smooth_sigma,
                include_absorption=include_absorption,
            )
        if index == 1 or index == len(sparcl_ids) or index % 25 == 0:
            print(f"  Processed {index}/{len(sparcl_ids)} spectra (latest SPARCL ID: {sparcl_id})")

    sentinel_path.touch()
    print(f"  Saved {len(sparcl_ids)} spectra to {output_dir}")


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Query DESI/SPARCL spectra around a sky position, download all matching "
            "spectra, and optionally save diagnostic plots."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --output clstr01
python desi_download_spectra.py --csv targets.csv
python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --no-plots
""",
    )
    parser.add_argument("--ra", type=float, help="Right ascension in decimal degrees.")
    parser.add_argument("--dec", type=dec_degrees, help="Declination in decimal degrees.")
    parser.add_argument("--radius", type=positive_float, help="Cone search radius in degrees.")
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for a single target (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help=(
            "Run the script in batch mode using a CSV file instead of a single --ra/--dec/--radius target. "
            "The CSV must contain columns: name, ra, dec, radius. An optional output column is also supported."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=DEFAULT_LIMIT,
        help=f"Maximum number of spectra returned per target (default: {DEFAULT_LIMIT}).",
    )
    parser.add_argument(
        "--data-release",
        default=DEFAULT_DATA_RELEASE,
        help=f"SPARCL data release filter (default: {DEFAULT_DATA_RELEASE}).",
    )
    parser.add_argument(
        "--table",
        default="sparcl.main",
        help="SPARCL table name to query (default: sparcl.main).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Reprocess a target even if download_complete.txt already exists.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip PNG spectrum plots and save only .npz files plus the object catalog.",
    )
    parser.add_argument(
        "--show-model",
        action="store_true",
        help="Overlay the SPARCL model spectrum on each plot when available.",
    )
    parser.add_argument(
        "--no-smooth",
        action="store_true",
        help="Disable the smoothed spectrum overlay in plots.",
    )
    parser.add_argument(
        "--smooth-sigma",
        type=positive_float,
        default=DEFAULT_SMOOTH_SIGMA,
        help=f"Gaussian smoothing sigma in pixels for plots (default: {DEFAULT_SMOOTH_SIGMA}).",
    )
    parser.add_argument(
        "--no-absorption-lines",
        action="store_true",
        help="Do not overlay common absorption lines on plots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.limit <= 0:
        raise SystemExit("--limit must be greater than 0")

    save_plots = not args.no_plots
    show_smooth = not args.no_smooth
    include_absorption = not args.no_absorption_lines

    if args.csv:
        csv_path = args.csv.expanduser().resolve()
        try:
            targets = read_targets_csv(csv_path)
        except (FileNotFoundError, ValueError) as exc:
            raise SystemExit(str(exc)) from exc

        client = create_sparcl_client()
        total_targets = len(targets)
        for target_index, row in enumerate(targets, start=1):
            output_name = row.get("output") or row["name"]
            process_target(
                client=client,
                ra_deg=float(row["ra"]),
                dec_deg=float(row["dec"]),
                radius_deg=float(row["radius"]),
                output_dir=Path(str(output_name)),
                target_label=str(row["name"]),
                target_index=target_index,
                total_targets=total_targets,
                table=args.table,
                data_release=args.data_release,
                limit=args.limit,
                overwrite=args.overwrite,
                save_plots=save_plots,
                show_model=args.show_model,
                show_smooth=show_smooth,
                smooth_sigma=args.smooth_sigma,
                include_absorption=include_absorption,
            )
        return

    if args.ra is None or args.dec is None or args.radius is None:
        raise SystemExit("--ra, --dec, and --radius are required unless --csv is provided.")

    client = create_sparcl_client()
    process_target(
        client=client,
        ra_deg=args.ra,
        dec_deg=args.dec,
        radius_deg=args.radius,
        output_dir=Path(args.output),
        target_label=args.output,
        target_index=1,
        total_targets=1,
        table=args.table,
        data_release=args.data_release,
        limit=args.limit,
        overwrite=args.overwrite,
        save_plots=save_plots,
        show_model=args.show_model,
        show_smooth=show_smooth,
        smooth_sigma=args.smooth_sigma,
        include_absorption=include_absorption,
    )


if __name__ == "__main__":
    main()
