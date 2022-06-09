# Astro Scripts
Various astronomy scripts

1. Download DECalS Image Script `decal_image_download.py`

    Downloading optical image from DeCALS survers ([https://www.legacysurvey.org/decamls/](https://www.legacysurvey.org/decamls/))

    `python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>`

    - \<redshift\>: set value to `0` when you do not know the redshift
    - \<fraction_size\>: 1 = 2048px, 2 = 2048/2, 4 = 2048/4; (Increase to zoom in and reduce the file size)

    Example:  
    Download jpeg: `python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4`
    
    Download FITS file: `python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits`
    

