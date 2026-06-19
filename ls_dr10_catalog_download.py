#!/usr/bin/env python3

import numpy as np
import argparse
from astroquery.utils.tap.core import TapPlus


def build_query(ra_min, ra_max, dec_min, dec_max):
    return f"""
    SELECT
        ra, dec,

        dered_flux_g, dered_flux_r, dered_flux_i, dered_flux_z,
        dered_flux_w1, dered_flux_w2,

        flux_ivar_g, flux_ivar_r, flux_ivar_i, flux_ivar_z,
        flux_ivar_w1, flux_ivar_w2,

        mw_transmission_g, mw_transmission_r, mw_transmission_i,
        mw_transmission_z, mw_transmission_w1, mw_transmission_w2,

        type, maskbits, shape_r,
        fracflux_g, fracflux_r, fracflux_i, fracflux_z

    FROM ls_dr10.tractor
    WHERE ra BETWEEN {ra_min} AND {ra_max}
      AND dec BETWEEN {dec_min} AND {dec_max}
    """


def main():
    parser = argparse.ArgumentParser(description="Download LS DR10 SED data")

    parser.add_argument("--name", type=str, default="cluster",
                        help="Object/field name for output file")

    parser.add_argument("--ra", type=float, required=True,
                        help="Central RA in degrees")

    parser.add_argument("--dec", type=float, required=True,
                        help="Central Dec in degrees")

    parser.add_argument("--radius", type=float, default=0.0166667,
                        help="Search radius in degrees")

    args = parser.parse_args()

    ra0 = args.ra
    dec0 = args.dec
    radius = args.radius

    # Compute bounding box
    dec_min = dec0 - radius
    dec_max = dec0 + radius

    ra_delta = radius / np.cos(np.deg2rad(dec0))
    ra_min = ra0 - ra_delta
    ra_max = ra0 + ra_delta

    # TAP connection
    tap = TapPlus(url="https://datalab.noirlab.edu/tap")

    query = build_query(ra_min, ra_max, dec_min, dec_max)

    print("Submitting query...")
    job = tap.launch_job(query)
    tbl = job.get_results()

    df = tbl.to_pandas()

    # Clean filename
    filename = f"ls_dr10_{args.name}_ra{ra0:.5f}_dec{dec0:.5f}_r{radius:.5f}.csv"

    df.to_csv(filename, index=False)

    print(f"Saved file: {filename}")
    print(f"Number of sources: {len(df)}")


if __name__ == "__main__":
    main()