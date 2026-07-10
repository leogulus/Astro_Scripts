#!/usr/bin/env python3
"""
download_desi_spectra.py

Download and plot DESI/SPARCL spectra within a cone centered on a given
(RA, Dec) position.

Usage
-----
Basic usage:
    python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02

Specify an output directory:
    python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --output clstr01

Arguments
---------
--ra
    Right Ascension of the search center in decimal degrees.

--dec
    Declination of the search center in decimal degrees.

--radius
    Cone search radius in decimal degrees.

--output
    Output directory (default: desi_output).

Outputs
-------
The script creates the output directory and saves:

    object_catalog.csv
        Catalog of all objects returned by the cone search.

    spectra.pkl
        Pickled SPARCL/specutils object returned by `client.retrieve()`.

    spectrum_<sparcl_id>.png
        One annotated spectrum plot for each downloaded object.

Dependencies
------------
Required Python packages:

    numpy
    matplotlib
    astropy
    sparcl
    dl

Example
-------
$ python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --output clstr01
"""

import os
import argparse
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path

from astropy.convolution import convolve, Gaussian1DKernel
from dl import queryClient as qc
from sparcl.client import SparclClient

client = SparclClient()

def cone_query(
    ra0: float,
    dec0: float,
    radius: float,
    table: str = "sparcl.main",
    extra_where: str = "data_release='DESI-DR1'",
    limit: int = 1000,
):
    """Perform a cone search in the SPARCL catalog."""

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

    cosdec = max(np.cos(np.radians(dec0)), 1e-6)
    ra_width = radius / cosdec
    ra_min = ra0 - ra_width
    ra_max = ra0 + ra_width
    dec_min = dec0 - radius
    dec_max = dec0 + radius

    sql = f"""
    SELECT {", ".join(cols)}
    FROM {table}
    WHERE
        ra BETWEEN {ra_min} AND {ra_max}
        AND dec BETWEEN {dec_min} AND {dec_max}
        AND acos(
            sin(radians(dec))*sin(radians({dec0})) +
            cos(radians(dec))*cos(radians({dec0}))*
            cos(radians(ra-{ra0}))
        ) <= radians({radius})
        AND {extra_where}
    LIMIT {limit}
    """

    return qc.query(sql=sql, fmt="pandas")


def get_default_lines(include_absorption=False):
    """Return default spectral lines for plotting."""

    lines = {
        "[O II]": 3727.0,
        "Hβ": 4861.33,
        "[O III]": 5006.84,
        "Hα": 6562.80,
        "[N II]": 6583.45,
        "[S II]": 6716.44,
    }

    if include_absorption:
        lines.update(
            {
                "Ca K": 3933.66,
                "Ca H": 3968.47,
                "Hδ": 4101.74,
                "Hγ": 4340.47,
                "Mg b": 5175.27,
                "Na D": 5895.92,
            }
        )

    return lines

def meta_value(meta, key, idx):
    value = meta.get(key)
    if value is None:
        return None
    try:
        return value[idx]
    except (TypeError, KeyError, IndexError):
        return value

def save_spectrum(results, idx, output_dir):
    """Save a single DESI spectrum to a compressed NumPy file."""

    sid = results.meta["sparcl_id"][idx]

    def get_meta(key):
        value = results.meta.get(key)
        if value is None:
            return None
        return np.asarray(value[idx])

    outfile = os.path.join(output_dir, f"spectrum_{sid}.npz")

    np.savez(
        outfile,
        wavelength=np.asarray(results.spectral_axis),
        flux=np.asarray(results.flux[idx]),
        model=get_meta("model"),
        ivar=get_meta("ivar"),
        mask=get_meta("mask"),
        wave_sigma=get_meta("wave_sigma"),
        sparcl_id=sid,
        redshift = meta_value(results.meta, "redshift", idx),
        specid = meta_value(results.meta, "specid", idx),
        spectype = meta_value(results.meta, "spectype", idx),
        ra = meta_value(results.meta, "ra", idx),
        dec = meta_value(results.meta, "dec", idx),
    )

    return outfile

def plot_spectrum(
    record,
    idx,
    output_dir,
    show_model=False,
    show_smooth=True,
    smooth_sigma=5,
    include_absorption=True,
):
    sid = record.meta["sparcl_id"][idx]
    outfile = os.path.join(output_dir, f"spectrum_{sid}.png")

    wave = np.asarray(record.spectral_axis)
    flux = np.asarray(record.flux[idx]) if np.ndim(record.flux) != 1 else record.flux
    wave_unit = getattr(record.spectral_axis, "unit", "")
    flux_unit = getattr(record.flux, "unit", "")

    meta = {k: meta_value(record.meta, k, idx) for k in record.meta}
    z = meta.get("redshift", 0.0)

    plt.figure(figsize=(10, 5))
    plt.plot(wave, flux, alpha=0.3, lw=0.8, label="Observed")

    if show_smooth:
        kernel = Gaussian1DKernel(smooth_sigma)
        smooth = convolve(flux, kernel, boundary="extend")
        plt.plot(wave, smooth, lw=1.2, label="Smoothed")

    if show_model and "model" in record.meta:
        model = np.asarray(record.meta["model"][idx])
        if model.ndim > 1:
            model = model[0]
        plt.plot(wave, model, lw=1.5, label="Model")

    ymax = plt.ylim()[1]

    for name, lam in get_default_lines(include_absorption).items():
        obs = lam * (1 + z)
        plt.axvline(obs, linestyle="--", alpha=0.3)
        plt.text(
            obs,
            ymax * 0.95,
            name,
            rotation=90,
            fontsize=8,
            ha="right",
            va="top",
        )

    wave_label = f"Wavelength [{wave_unit}]" if wave_unit else "Wavelength"
    flux_label = f"Flux [{flux_unit}]" if flux_unit else "Flux"
    plt.xlabel(wave_label)
    plt.ylabel(flux_label)

    plt.title(
        f"SPARCL ID = {meta['sparcl_id']}\n"
        f"z = {z:.4f}, "
        f"RA = {meta['ra']:.5f}, "
        f"Dec = {meta['dec']:.5f}"
    )

    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()

    return outfile

def get_meta_array(meta, key, index):
    value = meta.get(key)
    if value is None:
        return None
    return np.asarray(value[index])

def process_target(ra,dec,radius,output_dir,target_index=1,n_targets=1):
    output_dir_p = Path(output_dir)
    sentinel = output_dir_p / "download_complete.txt"
    
    if output_dir_p.exists() and sentinel.exists():
        print(f"Skipping {output_dir_p}: already exists.")
        return
    
    output_dir_p.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Target {target_index}/{n_targets}: "
        f"{output_dir} ===")
    print(f"RA={ra:.6f}, Dec={dec:.6f}, Radius={radius:.4f}")

    found = cone_query(ra, dec, radius)

    if len(found) == 0:
        print("  No objects found.")
        return

    catalog_file = os.path.join(output_dir, "object_catalog.csv")
    found.to_csv(catalog_file, index=False)

    include = [
        "sparcl_id",
        "specid",
        "data_release",
        "redshift",
        "flux",
        "wavelength",
        "model",
        "ivar",
        "mask",
        "spectype",
        "ra",
        "dec",
        "wave_sigma",
    ]

    ids = list(found["sparcl_id"])

    results = client.retrieve(
        uuid_list=ids,
        include=include,
        fmt="specutils",
    )

    for i in tqdm(range(len(ids)), desc=output_dir_p.name):
        save_spectrum(results, i, output_dir_p)
        plot_spectrum(results, i, output_dir_p)

    sentinel.touch()
    print(f"  Saved {len(ids)} spectra.")
    

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Query DESI/SPARCL spectra around a sky position, "
            "download all matching spectra, and save plots."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --output clstr01
python desi_download_spectra.py --csv file.csv
""",
    )

    parser.add_argument(
        "--ra",
        type=float,
        help="Right Ascension in decimal degrees.",
    )

    parser.add_argument(
        "--dec",
        type=float,
        help="Declination in decimal degrees.",
    )

    parser.add_argument(
        "--radius",
        type=float,
        help="Cone search radius in degrees.",
    )

    parser.add_argument(
        "--output",
        default="desi_output",
        help="Output directory (default: desi_output).",
    )

    parser.add_argument(
        "--csv",
        type=str,
        help="CSV file containing columns: name, ra, dec, radius.",
    )

    args = parser.parse_args()
    
    if args.csv:
        targets = pd.read_csv(args.csv)

        required = {"name", "ra", "dec", "radius"}
        missing = required - set(targets.columns)
        if missing:
            raise ValueError(
                f"CSV file is missing required columns: {sorted(missing)}"
            )

        n_targets = len(targets)

        for i, (_, row) in enumerate(targets.iterrows(), start=1):
            process_target(
                ra=float(row["ra"]),
                dec=float(row["dec"]),
                radius=float(row["radius"]),
                output_dir=str(row["name"]),
                target_index=i,
                n_targets=n_targets,
            )
    else:
        if args.ra is None or args.dec is None or args.radius is None:
            parser.error("--ra, --dec, and --radius are required unless --csv is provided.")
        process_target(
            ra=args.ra,
            dec=args.dec,
            radius=args.radius,
            output_dir=args.output,
        )

if __name__ == "__main__":
    main()