import urllib.request, os
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

import argparse

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def _fetch(outfile, ra, dec, size, pxscale=0.262):
    """
    Pixel scale for DeCal: 0.262"/px
    """
    url=f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&size={size}&"+\
    f"layer=ls-dr9&pixscale={pxscale}&bands=grz"
    print(f"Downloading: {url}")
    urllib.request.urlretrieve(url,outfile)
    
def fetch_image_decal(ra, dec, size):
    """Return the data array for the image of object type"""
    filename = 'tmp.jpg'
    _fetch(filename, ra, dec, size=size)

def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.tick_params(labelbottom=False, labelleft=False)
        ax.axis('off')

def decal_image(indx, ra, dec, name, redshift=np.nan, fraction_size=1):
    """
    indx: index number (can be number)
    redshift: used to calculate physical size on the scale
    fraction size: 1 = 2048, 2 = 2048/2, 4=2048/4 (to change the zoom of the image)
    """
    size=int(2048/fraction_size); fraction0=fraction_size/2
    filename = f'img_ix{indx:03}_annoted_{name}.jpg'
    try:
        fetch_image_decal(ra, dec, size=size)
        img = mpimg.imread('tmp.jpg')
        fig=plt.figure(figsize=(10,10))
        imgplot = plt.imshow(img)
        plt.scatter(size/2,size/2,marker='x', color='white', lw=0.5)

        x=50/fraction0; y=140/fraction0
        plt.annotate("1 arcmin",(x+40,y-10/fraction0),color='white',fontsize=13)
        if np.isnan(redshift):
            plt.annotate(f"{name}",(x,y-50/fraction0),color='white',fontsize=13)
            dsize="nan"
        else:
            plt.annotate(f"{name} z={redshift:.2f}",(x,y-50/fraction0),color='white',fontsize=13)
            dsize=(cosmo.angular_diameter_distance(redshift)*(1*u.arcmin).to(u.radian).value).to(u.kpc)
        
        plt.annotate(f"{dsize:.0f}",(x+40,y+17/fraction0),color='white',fontsize=13)
        plt.annotate("", xy=(x, y), xytext=(x+229, y), arrowprops=dict(arrowstyle="-", color='white'))

        format_axes(fig)
        plt.tight_layout()
        plt.savefig(filename, dpi=120)
        plt.close()
        os.remove('tmp.jpg')
    
    except: 
        print('No Data')
        pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('index', help="Index number")
    parser.add_argument('ra', help='ra')
    parser.add_argument('dec', help='dec')
    parser.add_argument('name', help='object name')
    parser.add_argument('redshift', help='redshift', default=np.nan)
    parser.add_argument('fraction_size', help='fraction_size', default=1)

    args = parser.parse_args()
    decal_image(int(args.index), float(args.ra), float(args.dec), args.name, redshift=float(args.redshift), fraction_size=float(args.fraction_size))


if __name__ == '__main__':
    main()

