# Astro Scripts
Various astronomy scripts

1. Download DECalS Image Script: `decal_image_download.py`

    Downloading optical image from DeCALS survers ([https://www.legacysurvey.org/decamls/](https://www.legacysurvey.org/decamls/))

    `python decal_image_download.py <index> <ra> <dec> <name> <redshift> <fraction_size>`

    - \<redshift\>: set value to `0` when you do not know the redshift
    - \<fraction_size\>: 1 = 2048px, 2 = 2048/2, 4 = 2048/4; (Increase to zoom in and reduce the file size)

    Example:  
    Download jpeg: `python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4`
    
    Download FITS file: `python decal_image_download.py 1 57.8554 -15.4054 test 0.2 4 --fits`


    1.1 Fixed the order of FITS images in order to open with ds9
    
    Requirements: Need to install `aplpy` for this code to work. 
    
    `python decal_image_download_fits_fixed.py <outfilename_from_first_step.fits>`
    
2. Download SDSS Catalog Data: `download_SDSS_catalog_withSciServer.ipynb`

    - This code can only be run on SciServer https://apps.sciserver.org/compute/, and cannot be run on your own machine.   
    - To get access to the SciServer, You will have to create an account at https://apps.sciserver.org/login-portal/.  
    - Please follow the installation steps in the Example-Notebooks to be able to run your own codes https://github.com/sciserver/Example-Notebooks#installation.
    - Once you create a persistent folder in `/Storage/<username>/persistent/`, upload this code to that folder. The script should work in this environment.
    Additional Information:
        - One of the main components of the code is using SQL to query the data we want. Here is a quick tutorial on how to use SQL: <https://www.w3schools.com/sql/>. 
        - Description of each column in each SDSS Table can be found here <https://skyserver.sdss.org/dr14/en/help/docs/tabledesc.aspx>.
