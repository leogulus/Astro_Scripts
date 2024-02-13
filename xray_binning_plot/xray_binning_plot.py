import numpy as np
import matplotlib.pyplot as plt
import os, sys, argparse, ast
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

"""
input:
- annuli.txt
- spt0417_xspec_kt_err.csv: run different binning with error command to get a good error estimate
- spt0417_xspec_kt_full.csv: one run with a full annuli to match with <annuli.txt> file
- redshift of the cluster (for scaling)

Call: python xray_binning_plot.py 'spt0417_xspec_kt_err.csv' 'spt0417_xspec_kt_full.csv' -z 0.58 -e
"""

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='spectra filename with error commands ran already', type=str)
    parser.add_argument('filename_full', help='spectra filename for all the binning without any cut (spt0417_xspec_kt_full.tsv)', type=str)
    parser.add_argument('-z','--redshift',help='redshift', default=None, type=float)
    parser.add_argument('-e','--error',help='already run the error command on XSPEC', action='store_true')
    args = parser.parse_args()
    return args

def read_parameter_file(filename):
    with open(filename) as f:
        lines = f.readlines()

    all_data={}
    for line in lines:
        line = line.strip()
        n=line.split()
        if len(n) == 1:
            data_i={}
            all_data[line]=data_i
        else:
            if 'err' in filename:
                error_num = ast.literal_eval(n[3])
                data_i[int(n[0])]=float(n[1])-error_num[0], np.abs(error_num[0]), error_num[1]
            else:
                try:
                    data_i[int(n[0])]=(float(n[5]), float(n[7]))
                except: 
                    data_i[int(n[0])]=(float(n[5]), n[7])
    return all_data

def match_annuli_fn(filename, all_data_full):
    with open(filename) as f:
        lines = f.readlines()

    all_annuli=[]
    for line in lines:
        line = line.strip()
        n=line.split()
        all_annuli.append(float(n[1]))

    match_annuli = dict(zip(all_data_full[list(all_data_full)[0]].keys(), all_annuli))
    return match_annuli


def create_new_binning(Alist, Blist):
    result_dict = {} 
    for i in range(len(Alist)):
        start_loc = Alist[i]
        if i == len(Alist)-1:
            end_loc = Blist[-1]
        else:
            end_loc = Alist[i+1]
        result_dict[Alist[i]] = ((end_loc+start_loc)/2, (end_loc+start_loc)/2-start_loc, end_loc-(end_loc+start_loc)/2)
    return [[i[0] for i in result_dict.values()], [i[1] for i in result_dict.values()], [i[2] for i in result_dict.values()]]

def main():
    args = parseArguments()
    
    all_data =read_parameter_file(args.filename)
    all_data_full = read_parameter_file(args.filename_full)

    scale_redshift=cosmo.kpc_proper_per_arcmin(args.redshift).to(u.kpc/u.arcsec).value #kpc/arcsec

    match_annuli = match_annuli_fn('annuli.txt', all_data_full)
    Blist = np.array([match_annuli[k]*0.492*scale_redshift for k in all_data_full[list(all_data_full)[0]].keys()])

    plt.figure(figsize=(4*2,8))
    for indx, bin_i in enumerate(all_data.keys()):
        plt.subplot(round(len(all_data.keys())/2)+1,2,indx+1)
        di = all_data[bin_i]
        
        Alist = np.array([match_annuli[k]*0.492*scale_redshift for k in di.keys()])
        newbin = create_new_binning(Alist,Blist)
        newbin[1][0] = newbin[0][0] #make the innermost bin all the way to zero
        if not args.error:
        # for first pass (w/o err function)
            plt.errorbar(newbin[0], [i[0] for i in di.values()], \
                         xerr=[newbin[1], newbin[2]], yerr=[i[1] for i in di.values()], fmt='o', label=bin_i)
        else:    
        # w/ err function
            plt.errorbar(newbin[0], [i[0] for i in di.values()], \
                         xerr=[newbin[1], newbin[2]], yerr=[[i[1] for i in di.values()],[i[2] for i in di.values()]], fmt='o', label=bin_i)

        plt.ylim(2,10)
        plt.xlim(2,7e2)
        plt.xlabel('d [kpc]')
        plt.ylabel('cl.apec.kT [keV]')
        plt.xscale('log')
        plt.legend()
    plt.tight_layout()
    
    filename, _ = os.path.splitext(args.filename)
    plt.savefig(f"{filename}_plot.jpg")

if __name__ == "__main__":
    main()