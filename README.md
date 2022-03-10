# Astro Scripts
Various astronomy scripts

1. Download DECalS Image Script `decal_image_download.py`

    Downloading optical image from DeCALS survers ([https://www.legacysurvey.org/decamls/](https://www.legacysurvey.org/decamls/))

    `python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>`
    
    - \<redshift\> and \<fraction_size\> are optional
    - \<fraction_size\>: 1 = 2048px, 2 = 2048/2, 4=2048/4; (Increase to zoom in)
