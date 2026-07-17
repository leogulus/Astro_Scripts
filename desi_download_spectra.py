#!/usr/bin/env python3
"""Download DESI/SPARCL spectra for one target or a batch of targets."""

from __future__ import annotations

import argparse
import csv
import math
import time
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
SPARCL_SELECT_COLUMNS = [
    "sparcl_id",
    "specid",
    "targetid",
    "ra",
    "dec",
    "redshift",
    "redshift_err",
    "spectype",
    "data_release",
]
DESI_JOIN_SELECT_COLUMNS = [
    "z.targetid",
    "z.z AS zpix_z",
    "z.zerr AS zpix_zerr",
    "z.zwarn",
    "z.chi2",
    "z.deltachi2",
    "z.subtype",
    "z.survey",
    "z.program",
    "z.healpix",
    "z.coadd_exptime",
    "z.coadd_numexp",
    "z.coadd_numnight",
    "z.coadd_numtile",
    "z.main_primary",
    "z.zcat_primary",
    "z.mean_fiber_ra",
    "z.mean_fiber_dec",
    "p.ls_id",
    "p.flux_g",
    "p.flux_r",
    "p.flux_z",
    "p.flux_w1",
    "p.flux_w2",
    "p.flux_ivar_g",
    "p.flux_ivar_r",
    "p.flux_ivar_z",
    "p.flux_ivar_w1",
    "p.flux_ivar_w2",
    "p.fiberflux_g",
    "p.fiberflux_r",
    "p.fiberflux_z",
    "p.fibertotflux_g",
    "p.fibertotflux_r",
    "p.fibertotflux_z",
    "p.ebv",
    "p.maskbits",
    "p.morphtype",
    "p.shape_r",
    "p.parallax",
    "p.pmra",
    "p.pmdec",
    "p.gaia_phot_g_mean_mag",
    "p.gaia_phot_bp_mean_mag",
    "p.gaia_phot_rp_mean_mag",
    "t.priority_init",
    "t.numobs_init",
    "t.obsconditions",
    "t.photsys",
    "t.desi_target AS target_desi_target",
    "t.bgs_target AS target_bgs_target",
    "t.mws_target AS target_mws_target",
    "t.scnd_target AS target_scnd_target",
]
JOIN_CHUNK_SIZE = 500
DEFAULT_RETRY_ATTEMPTS = 3
DEFAULT_RETRY_DELAY = 3.0
RETRIEVE_CHUNK_SIZE = 100
CATALOG_COLUMN_PRESETS = {
    "default": [
        "sparcl_id",
        "specid",
        "targetid",
        "ra",
        "dec",
        "redshift",
        "redshift_err",
        "zpix_z",
        "zpix_zerr",
        "zwarn",
        "spectype",
        "subtype",
        "chi2",
        "deltachi2",
        "survey",
        "program",
        "healpix",
        "coadd_numexp",
        "coadd_numnight",
        "coadd_numtile",
        "ls_id",
        "flux_g",
        "flux_r",
        "flux_z",
        "flux_w1",
        "flux_w2",
        "maskbits",
        "morphtype",
        "shape_r",
        "data_release",
    ],
    "duplicates": [
        "sparcl_id",
        "specid",
        "targetid",
        "ls_id",
        "ra",
        "dec",
        "mean_fiber_ra",
        "mean_fiber_dec",
        "redshift",
        "redshift_err",
        "zpix_z",
        "zpix_zerr",
        "zwarn",
        "chi2",
        "deltachi2",
        "spectype",
        "subtype",
        "survey",
        "program",
        "healpix",
        "coadd_exptime",
        "coadd_numexp",
        "coadd_numnight",
        "coadd_numtile",
        "main_primary",
        "zcat_primary",
        "flux_g",
        "flux_r",
        "flux_z",
        "flux_w1",
        "flux_w2",
        "maskbits",
        "morphtype",
        "shape_r",
        "parallax",
        "pmra",
        "pmdec",
        "target_desi_target",
        "target_bgs_target",
        "target_mws_target",
        "target_scnd_target",
        "data_release",
    ],
    "full": [],
}


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
    cos_dec = max(np.cos(np.radians(dec_deg)), 1e-6)
    ra_width = radius_deg / cos_dec
    ra_min = ra_deg - ra_width
    ra_max = ra_deg + ra_width
    dec_min = max(-90.0, dec_deg - radius_deg)
    dec_max = min(90.0, dec_deg + radius_deg)
    ra_clause = build_ra_clause(ra_min, ra_max)
    angular_distance = (
        f"acos(sin(radians(dec))*sin(radians({dec_deg:.10f})) + "
        f"cos(radians(dec))*cos(radians({dec_deg:.10f}))*cos(radians(ra-{ra_deg:.10f})))"
    )

    return f"""
    SELECT {", ".join(SPARCL_SELECT_COLUMNS)}
    FROM {table}
    WHERE
        {ra_clause}
        AND dec BETWEEN {dec_min:.10f} AND {dec_max:.10f}
        AND {angular_distance} <= radians({radius_deg:.10f})
        AND {extra_where}
    ORDER BY {angular_distance} ASC
    LIMIT {limit}
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


def run_with_retries(operation_name: str, func, *, attempts: int, delay_seconds: float):
    """Retry a network-backed operation a few times before failing."""
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            return func()
        except Exception as exc:
            last_error = exc
            if attempt == attempts:
                break
            print(
                f"  {operation_name} failed on attempt {attempt}/{attempts}: {exc}\n"
                f"  Retrying in {delay_seconds:.0f} seconds..."
            )
            time.sleep(delay_seconds)
    raise last_error


def data_release_schema(data_release: str) -> str:
    """Map a DESI data-release label to the corresponding Data Lab schema."""
    return data_release.lower().replace("-", "_")


def chunked(values: list[int], size: int) -> list[list[int]]:
    """Split a list into fixed-size chunks."""
    return [values[index:index + size] for index in range(0, len(values), size)]


def desi_join_sql(schema: str, targetid_chunk: list[int]) -> str:
    """Build a DESI metadata join for a chunk of target IDs."""
    targetids = ", ".join(str(targetid) for targetid in targetid_chunk)
    return f"""
    SELECT
        {", ".join(DESI_JOIN_SELECT_COLUMNS)}
    FROM {schema}.zpix AS z
    LEFT JOIN {schema}.photometry AS p
        ON z.targetid = p.targetid
    LEFT JOIN {schema}.target AS t
        ON z.targetid = t.targetid
    WHERE z.targetid IN ({targetids})
    """


def fetch_joined_desi_catalog(targetids, data_release: str):
    """Fetch additional DESI metadata from Data Lab tables."""
    from dl import queryClient as qc
    import pandas as pd

    valid_targetids = sorted(
        {
            int(targetid)
            for targetid in targetids
            if targetid is not None and str(targetid).strip().lower() not in {"", "nan", "none"}
        }
    )
    if not valid_targetids:
        return pd.DataFrame()

    schema = data_release_schema(data_release)
    frames = []
    for targetid_chunk in chunked(valid_targetids, JOIN_CHUNK_SIZE):
        sql = desi_join_sql(schema, targetid_chunk)
        frames.append(
            run_with_retries(
                "DESI metadata join query",
                lambda sql=sql: qc.query(sql=sql, fmt="pandas"),
                attempts=DEFAULT_RETRY_ATTEMPTS,
                delay_seconds=DEFAULT_RETRY_DELAY,
            )
        )

    joined = pd.concat(frames, ignore_index=True)
    return joined.drop_duplicates(subset=["targetid"], keep="first")


def build_final_catalog(found, data_release: str):
    """Merge the SPARCL cone-search result with richer DESI metadata."""
    joined = fetch_joined_desi_catalog(found["targetid"].tolist(), data_release=data_release)
    if len(joined) == 0:
        return found
    return found.merge(joined, on="targetid", how="left")


def deduplicate_by_targetid(found):
    """Keep one row per targetid, preferring the best-fit DESI solution."""
    if "targetid" not in found.columns:
        return found

    deduped = found.copy()
    if "chi2" in deduped.columns:
        deduped["chi2"] = deduped["chi2"].fillna(np.inf)
    if "redshift_err" in deduped.columns:
        deduped["redshift_err"] = deduped["redshift_err"].fillna(np.inf)
    if "deltachi2" in deduped.columns:
        deduped["deltachi2"] = deduped["deltachi2"].fillna(-np.inf)

    deduped = deduped.sort_values(
        by=["targetid", "redshift_err", "chi2", "deltachi2"],
        ascending=[True, True, True, False],
        na_position="last",
    )
    return deduped.drop_duplicates(subset=["targetid"], keep="first").reset_index(drop=True)


def select_catalog_columns(catalog, profile: str):
    """Apply a named output-column preset to the merged catalog."""
    if profile == "full":
        return catalog

    requested = CATALOG_COLUMN_PRESETS[profile]
    available = [column for column in requested if column in catalog.columns]
    return catalog.loc[:, available]


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
    if hasattr(found, "to_dict"):
        rows = found.to_dict(orient="records")
    else:
        rows = found

    for row in rows:
        sparcl_id = str(row["sparcl_id"])
        lookup[sparcl_id] = row
    return lookup


def read_catalog_csv(catalog_path: Path) -> list[dict[str, object]]:
    """Read a saved object_catalog.csv into a list of row dictionaries."""
    if not catalog_path.exists():
        raise FileNotFoundError(f"Catalog file not found: {catalog_path}")

    with catalog_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Catalog file is empty or missing a header row: {catalog_path}")
        return list(reader)


def save_spectrum(results, index: int, output_dir: Path) -> Path:
    """Save one spectrum as a compressed NumPy archive."""
    sparcl_id = results.meta["sparcl_id"][index]

    def array_meta(key: str):
        value = results.meta.get(key)
        if value is None:
            return None
        return np.asarray(value[index])

    output_path = output_dir / f"spectrum_{sparcl_id}.npz"
    payload = {
        "wavelength": np.asarray(results.spectral_axis),
        "flux": np.asarray(results.flux[index]) if np.ndim(results.flux) != 1 else np.asarray(results.flux),
        "sparcl_id": str(sparcl_id),
        "redshift": safe_float(metadata_value(results.meta, "redshift", index), default=float("nan")),
        "specid": str(metadata_value(results.meta, "specid", index) or ""),
        "spectype": str(metadata_value(results.meta, "spectype", index) or ""),
        "ra": safe_float(metadata_value(results.meta, "ra", index), default=float("nan")),
        "dec": safe_float(metadata_value(results.meta, "dec", index), default=float("nan")),
    }
    for key in ("model", "ivar", "mask", "wave_sigma"):
        value = array_meta(key)
        if value is not None:
            payload[key] = value
    np.savez(output_path, **payload)
    return output_path


def plot_spectrum_data(
    *,
    sparcl_id,
    wavelength,
    flux,
    model,
    wave_unit,
    flux_unit,
    metadata: dict[str, object],
    output_dir: Path,
    catalog_lookup: dict[str, dict[str, object]],
    show_model: bool,
    show_smooth: bool,
    smooth_sigma: float,
    include_absorption: bool,
) -> Path:
    """Plot one spectrum and save it to disk."""
    output_path = output_dir / f"spectrum_{sparcl_id}.png"

    catalog_row = catalog_lookup.get(str(sparcl_id), {})
    redshift = safe_float(metadata.get("redshift"), default=0.0)
    ra = safe_float(metadata.get("ra"), default=safe_float(catalog_row.get("ra"), default=float("nan")))
    dec = safe_float(metadata.get("dec"), default=safe_float(catalog_row.get("dec"), default=float("nan")))
    spectype = metadata.get("spectype") or catalog_row.get("spectype") or "unknown"
    subtype = catalog_row.get("subtype") or "unknown"
    survey = catalog_row.get("survey") or "unknown"
    program = catalog_row.get("program") or "unknown"
    morphtype = catalog_row.get("morphtype") or "unknown"

    plt.figure(figsize=(10, 5))
    plt.plot(wavelength, flux, alpha=0.3, lw=0.8, label="Observed")

    if show_smooth:
        kernel = Gaussian1DKernel(smooth_sigma)
        smoothed = convolve(flux, kernel, boundary="extend")
        plt.plot(wavelength, smoothed, lw=1.2, label="Smoothed")

    if show_model and model is not None:
        model = np.asarray(model)
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
        f"z = {redshift:.4f}, RA = {title_ra}, Dec = {title_dec}\n"
        f"{spectype} | {subtype} | {survey} | {program} | {morphtype}"
    )

    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()
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
    metadata = {key: metadata_value(results.meta, key, index) for key in results.meta}
    model = None
    if "model" in results.meta:
        model = results.meta["model"][index]

    return plot_spectrum_data(
        sparcl_id=sparcl_id,
        wavelength=np.asarray(results.spectral_axis),
        flux=np.asarray(results.flux[index]) if np.ndim(results.flux) != 1 else np.asarray(results.flux),
        model=model,
        wave_unit=getattr(results.spectral_axis, "unit", ""),
        flux_unit=getattr(results.flux, "unit", ""),
        metadata=metadata,
        output_dir=output_dir,
        catalog_lookup=catalog_lookup,
        show_model=show_model,
        show_smooth=show_smooth,
        smooth_sigma=smooth_sigma,
        include_absorption=include_absorption,
    )


def replot_directory(
    output_dir: Path,
    *,
    show_model: bool,
    show_smooth: bool,
    smooth_sigma: float,
    include_absorption: bool,
) -> None:
    """Regenerate PNG plots from saved NPZ spectra and object_catalog.csv."""
    catalog_path = output_dir / CATALOG_FILENAME
    catalog_rows = read_catalog_csv(catalog_path)
    catalog_lookup = build_catalog_lookup(catalog_rows)
    spectrum_files = sorted(output_dir.glob("spectrum_*.npz"))

    if not spectrum_files:
        raise FileNotFoundError(f"No spectrum_*.npz files found in {output_dir}")

    print(f"\n=== Replotting {output_dir} ===")
    for index, spectrum_path in enumerate(spectrum_files, start=1):
        with np.load(spectrum_path, allow_pickle=False) as data:
            sparcl_id = str(data["sparcl_id"])
            metadata = {
                "sparcl_id": sparcl_id,
                "redshift": data["redshift"].item() if np.ndim(data["redshift"]) == 0 else data["redshift"],
                "specid": data["specid"].item() if np.ndim(data["specid"]) == 0 else data["specid"],
                "spectype": data["spectype"].item() if np.ndim(data["spectype"]) == 0 else data["spectype"],
                "ra": data["ra"].item() if np.ndim(data["ra"]) == 0 else data["ra"],
                "dec": data["dec"].item() if np.ndim(data["dec"]) == 0 else data["dec"],
            }
            model = None
            if "model" in data.files:
                candidate = data["model"]
                if candidate.dtype != object or candidate.item() is not None:
                    model = candidate

            plot_spectrum_data(
                sparcl_id=sparcl_id,
                wavelength=np.asarray(data["wavelength"]),
                flux=np.asarray(data["flux"]),
                model=model,
                wave_unit="Angstrom",
                flux_unit="1e-17 erg / (Angstrom s cm2)",
                metadata=metadata,
                output_dir=output_dir,
                catalog_lookup=catalog_lookup,
                show_model=show_model,
                show_smooth=show_smooth,
                smooth_sigma=smooth_sigma,
                include_absorption=include_absorption,
            )
        if index == 1 or index == len(spectrum_files) or index % 25 == 0:
            print(f"  Replotted {index}/{len(spectrum_files)} spectra")


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


def parse_target_row(row: dict[str, str], row_number: int) -> tuple[str, float, float, float, Path]:
    """Validate one batch CSV row and return normalized target values."""
    try:
        name = row["name"].strip()
        if not name:
            raise ValueError("name is empty")
        ra = ra_degrees(row["ra"])
        dec = dec_degrees(row["dec"])
        radius = positive_float(row["radius"])
    except (KeyError, TypeError, ValueError, argparse.ArgumentTypeError) as exc:
        raise ValueError(f"row {row_number}: invalid target data ({exc})") from exc

    output_name = (row.get("output") or name).strip()
    if not output_name:
        raise ValueError(f"row {row_number}: output is empty")
    return name, ra, dec, radius, Path(output_name)


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
    catalog_columns: str,
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
        found = run_with_retries(
            "Cone search",
            lambda: cone_query(
                ra_deg,
                dec_deg,
                radius_deg,
                table=table,
                data_release=data_release,
                limit=limit,
            ),
            attempts=DEFAULT_RETRY_ATTEMPTS,
            delay_seconds=DEFAULT_RETRY_DELAY,
        )
    except Exception as exc:
        raise SystemExit(f"Cone search failed for {target_label}: {exc}") from exc

    if len(found) == 0:
        print("  No objects found.")
        return

    try:
        final_catalog = build_final_catalog(found, data_release=data_release)
    except Exception as exc:
        raise SystemExit(f"DESI catalog join failed for {target_label}: {exc}") from exc

    original_count = len(final_catalog)
    final_catalog = deduplicate_by_targetid(final_catalog)
    removed_count = original_count - len(final_catalog)
    if removed_count > 0:
        print(f"  Deduplicated {removed_count} rows with repeated targetid values using fit quality.")

    output_catalog = select_catalog_columns(final_catalog, catalog_columns)

    catalog_path = output_dir / CATALOG_FILENAME
    output_catalog.to_csv(catalog_path, index=False)

    sparcl_ids = [str(sparcl_id) for sparcl_id in final_catalog["sparcl_id"]]
    catalog_lookup = build_catalog_lookup(final_catalog)
    received_ids: set[str] = set()
    processed = 0
    for request_chunk in chunked(sparcl_ids, RETRIEVE_CHUNK_SIZE):
        try:
            results = run_with_retries(
                "Spectrum retrieval",
                lambda request_chunk=request_chunk: client.retrieve(
                    uuid_list=request_chunk,
                    include=include_fields(show_model),
                    fmt="specutils",
                ),
                attempts=DEFAULT_RETRY_ATTEMPTS,
                delay_seconds=DEFAULT_RETRY_DELAY,
            )
        except Exception as exc:
            raise SystemExit(f"Spectrum retrieval failed for {target_label}: {exc}") from exc

        returned_ids = [str(value) for value in results.meta.get("sparcl_id", [])]
        for result_index, sparcl_id in enumerate(returned_ids):
            if sparcl_id not in request_chunk or sparcl_id in received_ids:
                print(f"  Ignoring unexpected or duplicate returned SPARCL ID: {sparcl_id}")
                continue
            save_spectrum(results, result_index, output_dir)
            if save_plots:
                plot_spectrum(
                    results,
                    result_index,
                    output_dir,
                    catalog_lookup,
                    show_model=show_model,
                    show_smooth=show_smooth,
                    smooth_sigma=smooth_sigma,
                    include_absorption=include_absorption,
                )
            received_ids.add(sparcl_id)
            processed += 1
            if processed == 1 or processed == len(sparcl_ids) or processed % 25 == 0:
                print(f"  Processed {processed}/{len(sparcl_ids)} spectra (latest SPARCL ID: {sparcl_id})")

    missing_ids = set(sparcl_ids) - received_ids
    if missing_ids:
        print(f"  Incomplete retrieval: {len(missing_ids)} requested spectra were not returned; retry this target.")
        return

    sentinel_path.touch()
    print(f"  Saved {processed} spectra to {output_dir}")


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
    parser.add_argument("--ra", type=ra_degrees, help="Right ascension in decimal degrees.")
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
        "--replot-dirs",
        nargs="+",
        type=Path,
        help=(
            "Regenerate spectrum PNG files from existing output directories that already contain "
            "object_catalog.csv and spectrum_*.npz files, without doing any new search or download."
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
        "--catalog-columns",
        choices=["default", "duplicates", "full"],
        default="default",
        help=(
            "Choose which columns are written to object_catalog.csv: "
            "'default' for normal use, 'duplicates' for diagnosing repeated rows, "
            "or 'full' for every fetched join column."
        ),
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

    if args.replot_dirs:
        for output_dir in args.replot_dirs:
            replot_directory(
                output_dir.expanduser().resolve(),
                show_model=args.show_model,
                show_smooth=show_smooth,
                smooth_sigma=args.smooth_sigma,
                include_absorption=include_absorption,
            )
        return

    if args.csv:
        csv_path = args.csv.expanduser().resolve()
        try:
            targets = read_targets_csv(csv_path)
        except (FileNotFoundError, ValueError) as exc:
            raise SystemExit(str(exc)) from exc

        client = create_sparcl_client()
        total_targets = len(targets)
        for target_index, row in enumerate(targets, start=1):
            try:
                name, ra, dec, radius, output_dir = parse_target_row(row, target_index + 1)
            except ValueError as exc:
                print(f"Skipping CSV target: {exc}")
                continue
            process_target(
                client=client,
                ra_deg=ra,
                dec_deg=dec,
                radius_deg=radius,
                output_dir=output_dir,
                target_label=name,
                target_index=target_index,
                total_targets=total_targets,
                table=args.table,
                data_release=args.data_release,
                limit=args.limit,
                catalog_columns=args.catalog_columns,
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
        catalog_columns=args.catalog_columns,
        overwrite=args.overwrite,
        save_plots=save_plots,
        show_model=args.show_model,
        show_smooth=show_smooth,
        smooth_sigma=args.smooth_sigma,
        include_absorption=include_absorption,
    )


if __name__ == "__main__":
    main()
