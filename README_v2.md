# Astro Scripts
A small collection of astronomy utilities for image retrieval, catalog queries, and photometry processing.
---
## 1. DECaLS Image Download
**Script:** `decal_image_download.py`
Downloads optical images from the DECaLS survey  
https://www.legacysurvey.org/decamls/
### Usage
```bash
python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>

Inputs

* ra: Right Ascension (deg)
* dec: Declination (deg)
* name: Output filename prefix
* redshift: Use 0 if unknown
* fraction_size:
    * 1 → 2048 px
    * 2 → 1024 px
    * 4 → 512 px
        Larger value means higher zoom and smaller output size

Examples

Download JPEG:

python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4

Download FITS:

python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits

⸻

FITS Fix for DS9

Fix FITS extension ordering for proper DS9 visualization.

Requirement: aplpy

python decal_image_download_fits_fixed.py <input_fits_file>

⸻

2. SDSS Catalog Query (SciServer)

Notebook: download_SDSS_catalog_withSciServer.ipynb

This notebook runs on SciServer only
https://apps.sciserver.org/compute/

Setup

* Create an account: https://apps.sciserver.org/login-portal/
* Follow setup instructions: https://github.com/sciserver/Example-Notebooks#installation
* Run inside: /Storage/<username>/persistent/

Notes

* Uses SQL queries to access SDSS data
* SQL reference: https://www.w3schools.com/sql/
* SDSS DR14 schema: https://skyserver.sdss.org/dr14/en/help/docs/tabledesc.aspx

⸻

3. Legacy Survey DR10 Photometry

Script: ls_dr10_sed_download.py

Queries the Legacy Survey DR10 Tractor catalog via NOIRLab TAP
https://datalab.noirlab.edu/tap

Usage

python ls_dr10_sed_download.py --ra <ra> --dec <dec> --name <name> --radius <radius>

Inputs

* ra: Right Ascension (deg)
* dec: Declination (deg)
* name: Output filename prefix
* radius: Search radius in degrees (default: 0.0166667 ≈ 1 arcmin)

Output

A CSV file containing:

* Sky coordinates (RA, Dec)
* Optical and WISE fluxes
* Flux inverse variance
* Galactic transmission values
* Source morphology and quality flags

⸻

Photometry and SED Calculations

For each band (g, r, i, z, w1, w2, …):

Galactic extinction correction

flux_corr = flux / mw_transmission

Flux uncertainty

sigma_flux = 1 / sqrt(flux_ivar)
sigma_corr = sigma_flux / mw_transmission

Magnitude

mag = 22.5 - 2.5 * log10(flux_corr)

Magnitude uncertainty

mag_err = (2.5 / ln(10)) * (sigma_flux / flux)