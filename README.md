# Astro Scripts

A collection of Python scripts and notebooks for retrieving astronomical imaging and catalog data from major surveys (DECaLS, SDSS, and LS DR10).

---

## Table of Contents
* [1. LS DR10 Photometry Downloader (`ls_dr10_catalog_download.py`)](#1-ls-dr10-photometry-downloader-ls_dr10_catalog_downloadpy)
* [2. DESI Spectra Downloader (`desi_download_spectra.py`)](#2-desi-spectra-downloader-desi_download_spectrapy)
* [3. DECaLS Image Downloader (`decal_image_download.py`)](#3-decals-image-downloader-decal_image_downloadpy)
* [4. SDSS Catalog Downloader (`download_SDSS_catalog_withSciServer.ipynb`)](#4-sdss-catalog-downloader-download_sdss_catalog_withsciserveripynb)

---

## 1. LS DR10 Photometry Downloader (`ls_dr10_catalog_download.py`)

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
* `--extra-columns`: Comma-separated list of additional `ls_dr10.tractor` columns to append to the default export.
* `--list-columns`: Query `TAP_SCHEMA.columns`, print the currently available `ls_dr10.tractor` columns, and exit.

### Example

Query sources around a target galaxy cluster:

```bash
python ls_dr10_catalog_download.py --ra 64.3950417 --dec -11.9110306 --name macs0417 --radius 0.02
```

Save the result to a specific directory:

```bash
python ls_dr10_catalog_download.py --ra 64.3950417 --dec -11.9110306 --name macs0417 --radius 0.02 --output-dir output_catalogs
```

Inspect which other columns can be downloaded:

```bash
python ls_dr10_catalog_download.py --list-columns
```

Append a few extra Tractor columns to the default CSV:

```bash
python ls_dr10_catalog_download.py --ra 64.3950417 --dec -11.9110306 --name macs0417 --radius 0.02 --extra-columns objid,brickid,ebv
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
* **Photometry:** Native observed fluxes (`flux_*`), observed magnitudes (`mag_*`), dereddened fluxes (`dered_flux_*`), and dereddened magnitudes (`dered_mag_*`)
* **Uncertainties / Quality:** Flux inverse variances (`flux_ivar_*`) and signal-to-noise ratios (`snr_*`)
* **Metadata:** Galactic transmission factors, morphology types, and quality flags

If `--extra-columns` is used, those columns are appended after the default set. The script validates requested names against the live `ls_dr10.tractor` schema before running the science query.

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

## 2. DESI Spectra Downloader (`desi_download_spectra.py`)

Queries the SPARCL catalog for DESI spectra around one sky position or a batch of positions, downloads the spectra, and saves both data products and diagnostic plots.

The script now uses a two-step metadata workflow:

1. Cone search in `sparcl.main` to find nearby spectra and collect `sparcl_id` / `targetid`
2. Join richer DESI metadata from `desi_dr1.zpix`, `desi_dr1.photometry`, and `desi_dr1.target`

The spectra themselves are still downloaded from SPARCL using `sparcl_id`.

### Requirements

The DESI script is heavier than the other scripts in this repository because it depends on two astronomy-specific client packages:

* `sparclclient`: provides `from sparcl.client import SparclClient`
* `astro-datalab`: provides `from dl import queryClient as qc`

These packages can be harder to install into an existing Python environment because they bring in a fairly specific scientific Python stack.

At the time of writing:

* `astro-datalab 2.24.0` declares `Requires-Python: >=3.9,<3.12`
* `astro-datalab 2.24.0` also pins versions of core packages such as `numpy`, `astropy`, `pandas`, `matplotlib`, `scipy`, `specutils`, `pyvo`, and `pycurl`
* `sparclclient 1.3.0` depends on `numpy`, `pandas`, `requests`, `specutils`, `spectres`, and `pyjwt`

Because of that, the safest approach is to install this script in its own clean conda environment instead of mixing it with unrelated packages.

### Recommended Conda Setup

Recommended starting point:

```bash
conda create -n desi python=3.11 pip
conda activate desi
python -m pip install --upgrade pip
python -m pip install astro-datalab==2.24.0 sparclclient==1.3.0
```

A tested working environment on this project used:

* `Python 3.11.15`
* `astro-datalab 2.24.0`
* `sparclclient 1.3.0`
* `numpy 1.26.3`
* `astropy 5.3.4`
* `pandas 2.1.4`
* `matplotlib 3.8.4`
* `specutils 1.13.0`

If you already have a complicated astronomy or machine-learning environment, it is strongly recommended not to install these DESI dependencies into that same environment.

### Single-Target Usage

```bash
python desi_download_spectra.py --ra <ra> --dec <dec> --radius <radius> --output <output_dir>
```

Example:

```bash
python desi_download_spectra.py --ra 64.3950417 --dec -11.9110306 --radius 0.02 --output macs0417
```

### Batch Usage With CSV

Use `--csv` to run the same workflow for many targets at once.

```bash
python desi_download_spectra.py --csv targets.csv
```

The included example file [targets.csv](/Users/taweewat/Documents/Softwares/Astro_Scripts/targets.csv) uses this format:

```csv
name,ra,dec,radius,output
macs0417,64.3950417,-11.9110306,0.02,macs0417
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
* `--catalog-columns default`: Write a compact everyday catalog.
* `--catalog-columns duplicates`: Write a catalog aimed at investigating repeated rows and target-level duplicates.
* `--catalog-columns full`: Write every fetched join column from the DESI tables.

### Duplicate Handling

If multiple rows share the same `targetid`, the script keeps only one row before downloading spectra.

The current preference order is:

* smaller `redshift_err`
* smaller `chi2`
* larger `deltachi2`

This helps reduce obvious repeated entries while keeping the better-fit DESI solution.

### Replot Existing Directories

If you already have output directories with `object_catalog.csv` and `spectrum_*.npz`, you can regenerate the PNG plots without doing the cone search, DESI joins, or SPARCL download again.

```bash
python desi_download_spectra.py --replot-dirs macs0417
```

This mode also works with plot options such as:

```bash
python desi_download_spectra.py --replot-dirs macs0417 --show-model --no-smooth
```

### Output

For each target, the script writes:

* `object_catalog.csv`: Final catalog after SPARCL cone search, DESI table joins, and optional column filtering
* `spectrum_<sparcl_id>.npz`: Saved spectrum arrays and metadata
* `spectrum_<sparcl_id>.png`: Plot of the spectrum, unless `--no-plots` is used
* `download_complete.txt`: Sentinel file marking the target as complete

---

## 3. DECaLS Image Downloader (`decal_image_download.py`)

Downloads optical images directly from the [DECaLS servers](https://www.legacysurvey.org/decamls/).

### Requirements

* `numpy`
* `matplotlib`
* `astropy`

### Usage

```bash
python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>
```

Both `<redshift>` and `<fraction_size>` are optional:

* default redshift: `0.2`
* default fraction size: `1`

### Arguments

* `<index>`: Identifier index or row number.
* `<ra>`: Right Ascension (in degrees).
* `<dec>`: Declination (in degrees).
* `<name>`: Output file name/prefix.
* `<redshift>`: Target redshift used for the physical scale label *(default: `0.2`)*. Set to `0` if unknown.
* `<fraction_size>`: Resolution scaling factor *(default: `1`)*.
* `--fits`: Download a FITS cube instead of an annotated JPEG.
* `--keep-raw-fits`: Keep the original downloaded FITS file alongside the reordered output.
* `--catalog-csv`: Path to `object_catalog.csv` from `desi_download_spectra.py`. When provided, the JPEG output marks those catalog objects on the image.
* `--brightest-csv`: Path to a CSV from `ls_dr10_catalog_download.py`. When provided, the JPEG output marks the brightest LS objects using `mag_i`.
* `--brightest-count`: Number of brightest LS objects to mark from `--brightest-csv` *(default: `5`)*.
* `--label-chars`: Number of leading `sparcl_id` characters to display beside each marker *(default: `3`)*.

`fraction_size` examples:
`1` = 2048 px
`2` = 2048 / 2 px
`4` = 2048 / 4 px *(Increase value to zoom in and reduce file size)*

### Examples

* **Download JPEG:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4
```

* **Download JPEG using the default redshift (`0.2`) and default image size:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417
```

* **Download FITS file:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4 --fits
```

* **Download FITS file and keep the original raw download too:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4 --fits --keep-raw-fits
```

* **Download a JPEG and overlay DESI catalog objects from an existing catalog CSV:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4 --catalog-csv macs0417/object_catalog.csv
```

* **Download a JPEG and overlay the 5 brightest LS DR10 sources from a photometry CSV:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4 --brightest-csv macs0417/ls_dr10_macs0417_ra64.39504_dec-11.91103_r0.02000.csv
```

* **Download a JPEG with both DESI and LS overlays at the same time:**

```bash
python decal_image_download.py 1 64.3950417 -11.9110306 macs0417 0.44 4 --catalog-csv macs0417/object_catalog.csv --brightest-csv macs0417/ls_dr10_macs0417_ra64.39504_dec-11.91103_r0.02000.csv
```

### Output

* JPEG mode writes an annotated image named like `img_ix00001_annoted_test.jpg`
* JPEG mode adds suffixes when overlays are used:
  `..._desi.jpg`, `..._lsbright.jpg`, or `..._desi_lsbright.jpg`
* FITS mode writes a reordered FITS cube named like `img_ix00001_annoted_test.fits`

The original FITS cube downloaded from the DECaLS server has the band order wrong for SAOImage DS9 RGB display. This script automatically reorders the downloaded cube before saving the final FITS output, so no second fix-up script is needed.

The JPEG annotation includes a real 1 arcmin scale bar computed from the DECaLS pixel scale (`0.262` arcsec/px), plus a physical-size label in kpc at the chosen redshift.

### DESI Marker Overlay

When `--catalog-csv` is used, the script reads that `object_catalog.csv` file and overlays the catalog objects on the JPEG cutout.

Each DESI marker is drawn as a hollow cyan circle and labeled with:

* the first `--label-chars` characters of `sparcl_id`
* the object redshift as `z=...` with 2 decimal places, when the `redshift` column is available in `object_catalog.csv`

The marker positions are computed using:

* the input `ra` and `dec` as the image center
* the DECaLS JPEG pixel scale of `0.262` arcsec/px

This mode does not use WCS; it uses sky-coordinate offsets from the image center and converts them directly into pixel offsets.

### LS Brightest Overlay

When `--brightest-csv` is used, the script reads the LS DR10 CSV, sorts by `mag_i`, and marks the brightest rows on the JPEG cutout.

Each LS marker is drawn as a hollow yellow circle and labeled with `mag_i` to 1 decimal place.

### DS9 Note

If you want to open the FITS cube in SAOImage DS9 as an RGB image, first go to `Frame -> RGB`, then choose `File -> Open As -> RGB Cube`.

---

## 4. SDSS Catalog Downloader (`download_SDSS_catalog_withSciServer.ipynb`)

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
