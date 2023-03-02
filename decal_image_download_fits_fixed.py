import os, argparse
from astropy.io import fits
import aplpy

def make_one_band_fits(filename,index,band):
    hdu = fits.open(filename)
    hdu[0].data=hdu[0].data[index]
    hdu.writeto(f"tmp_img_{band}.fits",overwrite=True)

def main():
    """
    Fixed the order of the image so that the RGB with ds9 makes sense
    output: 
        - outfilename: 
            - read more info here (https://aplpy.readthedocs.io/en/stable/api/aplpy.make_rgb_cube.html) 
    Example: 
        > python decal_image_download_fits_fixed.py img_ix001_annoted_SPT-CLJ0234-5831.fits
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('fitsname', help="fits filename")
    
    args = parser.parse_args()
    fitsname=args.fitsname
    outfilename=os.path.splitext(fitsname)[0]+"_fixed.fits"
    
    bands = ['g','r','z']
    for index, band in enumerate(bands):
        make_one_band_fits(fitsname,index,band)
        
    aplpy.make_rgb_cube(['tmp_img_z.fits','tmp_img_r.fits','tmp_img_g.fits'], outfilename, north=True)
    
    for band in bands:
        os.remove(f"tmp_img_{band}.fits")

if __name__ == '__main__':
    main()
