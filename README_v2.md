# Astro Scripts

A collection of Python scripts and notebooks for retrieving astronomical imaging and catalog data from major surveys (DECaLS, SDSS, and LS DR10).

---

## Table of Contents
* [1. DECaLS Image Downloader (`decal_image_download.py`)](#1-decals-image-downloader-decal_image_downloadpy)
* [2. SDSS Catalog Downloader (`download_SDSS_catalog_withSciServer.ipynb`)](#2-sdss-catalog-downloader-download_sdss_catalog_withsciserveripynb)
* [3. LS DR10 Photometry Downloader (`ls_dr10_catalog_download.py`)](#3-ls-dr10-photometry-downloader-ls_dr10_catalog_downloadpy)

---

## 1. DECaLS Image Downloader (`decal_image_download.py`)

Downloads optical images directly from the [DECaLS servers](https://www.legacysurvey.org/decamls/).

### Usage
```bash
python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size> [options]

```

### Arguments

* `<index>`: Identifier index or row number.
* `<ra>`: Right Ascension (in degrees).
* `<dec>`: Declination (in degrees).
* `<name>`: Output file name/prefix.
* `<redshift>`: Target redshift. Set to `0` if unknown.
* `<fraction_size>`: Resolution scaling factor.
* `1` = 2048 px
* `2` = $2048 / 2$ px
* `4` = $2048 / 4$ px *(Increase value to zoom in and reduce file size)*



### Examples

* **Download JPEG:**
```bash
python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4

```


* **Download FITS file:**
```bash
python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits

```



### 1.1 FITS Header Alignment (`decal_image_download_fits_fixed.py`)

Fixes the axis ordering of downloaded DECaLS FITS images so they can be parsed correctly by SAOImage DS9.

* **Requirements:** `aplpy`
* **Usage:**
```bash
python decal_image_download_fits_fixed.py <outfilename_from_first_step.fits>

```



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

Queries the Legacy Surveys DR10 Tractor catalog via the [NOIRLab TAP service](https://datalab.noirlab.edu/tap) to extract multi-band photometry (optical + WISE) within a user-defined sky region.

### Usage

```bash
python ls_dr10_catalog_download.py --ra <ra> --dec <dec> --name <name> --radius <radius>

```

### Arguments

* `--ra`: Right Ascension of the target center (in degrees).
* `--dec`: Declination of the target center (in degrees).
* `--name`: Label for the output CSV file.
* `--radius`: Search radius in degrees *(default: `0.0166667`, approximately 1 arcmin)*.

### Example

Query sources around a target galaxy cluster:

```bash
python ls_dr10_catalog_download.py --ra 150.11632 --dec 2.20583 --name clusterA --radius 0.02

```

### Output Columns

The script returns a CSV file containing:

* **Astrometry:** `ra`, `dec`
* **Photometry:** Optical and WISE fluxes (dereddened) and flux inverse variances (`flux_ivar`)
* **Metadata:** Galactic transmission factors, morphology types, and quality flags

### SED Photometry Calculations

For each filter band ($g, r, i, z, W1, W2, \dots$), the script processes the Tractor catalog inputs using the following conventions:

* **Dereddened Flux:**

$$f_{\text{dered}} = \frac{f}{\text{mw\_transmission}}$$


* **Flux Error:**

$$\sigma_{f} = \frac{1}{\sqrt{f_{\text{ivar}}}}$$


$$\sigma_{f, \text{dered}} = \frac{\sigma_{f}}{\text{mw\_transmission}}$$


* **AB Magnitude:**

$$m = 22.5 - 2.5 \log_{10}(f_{\text{dered}})$$


* **Magnitude Error:**

$$\sigma_{m} = \frac{2.5}{\ln(10)} \cdot \frac{1}{f \cdot \sqrt{f_{\text{ivar}}}}$$
```
