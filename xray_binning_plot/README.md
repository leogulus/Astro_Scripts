# X-ray Binning Plot

Input:
- `annuli.txt` (from *projct* folder)
- `spt0417_xspec_kt_err.csv`: run different binning with error command to get a good error estimate
- `spt0417_xspec_kt_full.csv`: one run with a full annuli to match with <annuli.txt> file
- redshift of the cluster (for arcsec to kpc scaling)

Call: `python xray_binning_plot.py 'spt0417_xspec_kt_err.csv' 'spt0417_xspec_kt_full.csv' -z 0.58 -e`

Output:
- One plot with all the binning 

![plot1](https://github.com/leogulus/Astro_Scripts/xray_binning_plot/spt0417_xspec_kt_err_plot.jpg?raw=true)
