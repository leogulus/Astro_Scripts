# Astro Scripts

A collection of Python scripts and notebooks for retrieving astronomical imaging and catalog data from major surveys (DECaLS, SDSS, and LS DR10).

---

## Table of Contents
* [1. DECaLS Image Downloader (`decal_image_download.py`)](#1-decals-image-downloader-decal_image_downloadpy)
* [2. SDSS Catalog Downloader (`download_SDSS_catalog_withSciServer.ipynb`)](#2-sdss-catalog-downloader-download_sdss_catalog_withsciserveripynb)
* [3. LS DR10 Photometry Downloader (`ls_dr10_catalog_download.py`)](#3-ls-dr10-photometry-downloader-ls_dr10_catalog_downloadpy)
* [4. DESI Spectra Downloader (`desi_download_spectra.py`)](#4-desi-spectra-downloader-desi_download_spectrapy)

---

## 1. DECaLS Image Downloader (`decal_image_download.py`)

Downloads optical images directly from the [DECaLS servers](https://www.legacysurvey.org/decamls/).

### Requirements

* `numpy`
* `matplotlib`
* `astropy`

### Usage
```bash
python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>
```

### Arguments

* `<index>`: Identifier index or row number.
* `<ra>`: Right Ascension (in degrees).
* `<dec>`: Declination (in degrees).
* `<name>`: Output file name/prefix.
* `<redshift>`: Target redshift. Set to `0` if unknown.
* `<fraction_size>`: Resolution scaling factor.
* `--fits`: Download a FITS cube instead of an annotated JPEG.
* `--keep-raw-fits`: Keep the original downloaded FITS file alongside the reordered output.

`fraction_size` examples:
`1` = 2048 px
`2` = 2048 / 2 px
`4` = 2048 / 4 px *(Increase value to zoom in and reduce file size)*

### Examples

* **Download JPEG:**
```bash
python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4
```

* **Download FITS file:**
```bash
python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits
```

* **Download FITS file and keep the original raw download too:**
```bash
python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits --keep-raw-fits
```

### Output

* JPEG mode writes an annotated image named like `img_ix00001_annoted_test.jpg`
* FITS mode writes a reordered FITS cube named like `img_ix00001_annoted_test.fits`

The original FITS cube downloaded from the DECaLS server has the band order wrong for SAOImage DS9 RGB display. This script automatically reorders the downloaded cube before saving the final FITS output, so no second fix-up script is needed.

### DS9 Note

If you want to open the FITS cube in SAOImage DS9 as an RGB image, first go to `Frame -> RGB`, then choose `File -> Open As -> RGB Cube`.

---

## 2. SDSS Catalog Downloader (`download_SDSS_catalog_withSciServer.ipynb`)

A Jupyter Notebook designed to query and download SDSS catalog data using the SciServer environment.

> [!IMPORTANT]
> This notebook relies on SciServer internal libraries and **cannot be run locally**. It must be executed within the SciServer Compute environment.

### Setup Instructions

1. Create an account at the [SciServer Login Portal](https://apps.sciserver.org/login-portal/).
2. Follow the setup steps in the [SciServer Example-Notebooks Installation Guide](https://github.com/sciserver/Example-Notebooks#installation).
3. Create a persistent directory under `/Storage/<username>/persistent/`.
4. Upload this notebook into that folder and run it within a SciServer Compute container.

### Useful Resources 

* **SQL Querying:** Data retrieval is handled via SQL. See this [W3Schools SQL Tutorial](https://www.w3schools.com/sql/) for a quick refresher.
* **SDSS Schema:** Column descriptions for all SDSS database tables are available at the [SkyServer DR14 Table Descriptions](https://skyserver.sdss.org/dr14/en/help/docs/tabledesc.aspx).

---

## 3. LS DR10 Photometry Downloader (`ls_dr10_catalog_download.py`)

Queries the Legacy Surveys DR10 Tractor catalog via the [NOIRLab TAP service](https://datalab.noirlab.edu/tap) to extract multi-band photometry (optical + WISE) within a user-defined sky box. 

More Information about LS DR10: https://datalab.noirlab.edu/data/legacy-surveys

### Requirements

* `astroquery`

### Search Method

This script intentionally performs a fast RA/Dec box search instead of a true cone search.
The returned table may include some sources near the corners of the box that fall outside a circular radius, which can be filtered manually afterward if needed.


### Usage

```bash
python ls_dr10_catalog_download.py --ra <ra> --dec <dec> --name <name> --radius <radius> --output-dir <output_dir>
```

### Arguments

* `--ra`: Right Ascension of the target center (in degrees).
* `--dec`: Declination of the target center (in degrees).
* `--name`: Label for the output CSV file.
* `--radius`: Half-width of the search box in degrees *(default: `0.0166667`, approximately 1 arcmin)*.
* `--output-dir`: Directory where the CSV file will be written *(default: current directory)*.

### Example

Query sources around a target galaxy cluster:

```bash
python ls_dr10_catalog_download.py --ra 150.11632 --dec 2.20583 --name clusterA --radius 0.02
```

Save the result to a specific directory:

```bash
python ls_dr10_catalog_download.py --ra 150.11632 --dec 2.20583 --name clusterA --radius 0.02 --output-dir output_catalogs
```

### Output

The script writes a CSV file named like:

```text
ls_dr10_<name>_ra<ra>_dec<dec>_r<radius>.csv
```

It also prints the RA/Dec bounds of the search box and the number of sources returned.

### Output Columns

The script returns a CSV file containing:

* **Astrometry:** `ra`, `dec`
* **Photometry:** Optical and WISE fluxes (dereddened) and flux inverse variances (`flux_ivar`)
* **Metadata:** Galactic transmission factors, morphology types, and quality flags

### SED Photometry Calculations

For each filter band (`g`, `r`, `i`, `z`, `w1`, `w2`, ...), the script converts Tractor catalog parameters using the following formulas:

```python
# 1. Use dereddened flux for SED
dered_flux = flux / mw_transmission

# 2. Flux error
sigma_flux = 1 / sqrt(flux_ivar)
sigma_dered_flux = sigma_flux / mw_transmission

# 3. AB Magnitude
mag = 22.5 - 2.5 * log10(dered_flux)

# 4. Magnitude error
mag_err = (2.5 / ln(10)) / (flux * sqrt(flux_ivar))
```

---

## 4. DESI Spectra Downloader (`desi_download_spectra.py`)

Queries the SPARCL catalog for DESI spectra around one sky position or a batch of positions, downloads the spectra, and saves both data products and diagnostic plots.

### Requirements

* `numpy`
* `matplotlib`
* `astropy`
* `sparcl`
* `dl`

### Single-Target Usage

```bash
python desi_download_spectra.py --ra <ra> --dec <dec> --radius <radius> --output <output_dir>
```

Example:

```bash
python desi_download_spectra.py --ra 140.1704 --dec 2.7832 --radius 0.02 --output clstr01
```

### Batch Usage With CSV

Use `--csv` to run the same workflow for many targets at once.

```bash
python desi_download_spectra.py --csv targets.csv
```

The included example file [targets.csv](/Users/taweewat/Documents/Softwares/Astro_Scripts/targets.csv) uses this format:

```csv
name,ra,dec,radius,output
cluster_a,140.1704,2.7832,0.02,clstr01
cluster_b,150.1163,2.2058,0.03,clstr02
```

Required CSV columns:

* `name`
* `ra`
* `dec`
* `radius`

Optional CSV column:

* `output`: Output directory name for that row. If omitted, the script uses `name`.

### Useful Options

* `--overwrite`: Re-run a target even if `download_complete.txt` already exists.
* `--no-plots`: Save only the catalog and `.npz` spectra, without PNG plots.
* `--show-model`: Overlay the SPARCL model spectrum when available.
* `--no-smooth`: Disable the smoothed spectrum overlay.
* `--smooth-sigma`: Change the Gaussian smoothing width for plots.
* `--no-absorption-lines`: Do not draw common absorption lines on the plot.

### Output

For each target, the script writes:

* `object_catalog.csv`: Catalog rows returned by the cone search
* `spectrum_<sparcl_id>.npz`: Saved spectrum arrays and metadata
* `spectrum_<sparcl_id>.png`: Plot of the spectrum, unless `--no-plots` is used
* `download_complete.txt`: Sentinel file marking the target as complete
