import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.optimize import curve_fit
import pickle

def gaussian( x, *p):
    c, mu, sigma = p
    res =   c * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    return res

def double_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2 = p
    res = A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))
    return res

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

def print_stats():
    print "\ngamma gaussian fit:\n  mu1 = "+str(coeff[1])+" ns\n  sigma1 = "+str(coeff[2]
            )+" ns"+'\n  mu2 = '+str(coeff[4])+'\n  sigma2 = '+str(coeff[5])    

save_dir = 'C:/Users/raweldon/Research/TUNL/beam_char_runs/10_17_run/analysis/0deg_tof/plots/'
# 11 MeV
# dists = ['180', '256', '363']
# 4 MeV
dists = ['179','276','369']

# get vaules by inspection with plot_tof_hist.py
# 11.325 MeV
#n_ranges=[[340.,342.],[323.5,326.],[299.8,302.]]
#g_ranges=[[400.,405.],[397.5,402.2],[394.,398.5]]
#n_p0s = [[1.0, 341., 1.0],[1.0, 324., 1.0],[1.0,300.,1.0]]
#g_p0s = [[1.0, 401., 0.1, 1.0, 403.5, 0.1],[1.0, 398.3, 0.1, 1.0, 401., 0.1],[1.0, 394.7, 0.1, 1.0, 397.3, 0.1]]
# 4.8 MeV
n_ranges=[[302.,307.],[269,275.],[240.5,246.]]
g_ranges=[[390.,410.],[398.,402.2],[394.,398.5]]
n_p0s = [[1.0, 302.5, 1.0],[1.0, 270., 1.0],[1.0,241.,1.0]]
g_p0s = [[1.0, 401., 0.1, 1.0, 403.5, 0.1],[1.0, 398.3, 0.1, 1.0, 400., 0.1],[1.0, 394.7, 0.1, 1.0, 397.3, 0.1]]

means_stds=[]

for index,dist in enumerate(dists):
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = get_range(tof, 200, 450)
    n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
    g_tof = get_range(tof, g_ranges[index][0], g_ranges[index][1])
    
    # build hists  
    tof_hist, bin_centers = build_hist(tof, 1e3)
    n_tof_hist, n_bin_centers = build_hist(n_tof, 500)
    g_tof_hist, g_bin_centers = build_hist(g_tof, 500)
    
    # neutrons
    p0 = n_p0s[index]
    coeff, var_matrix = curve_fit(gaussian, n_bin_centers, n_tof_hist, p0=p0)
    n_hist_fit = gaussian(n_bin_centers, *coeff)
    n_full_hist_fit = gaussian(bin_centers, *coeff)
    means_stds.append((coeff[1],coeff[2]))
    print "\nneutron gaussian fit:\n  mu = "+str(coeff[1])+" ns\n  sigma = "+str(coeff[2])

    # gammas
    p0 = g_p0s[index]
    coeff, var_matrix = curve_fit(double_gauss, g_bin_centers, g_tof_hist, p0=p0)
    g_hist_fit = double_gauss(g_bin_centers, *coeff)
    g_full_hist_fit = double_gauss(bin_centers, *coeff)
    means_stds.append((coeff[1],coeff[2],coeff[4],coeff[5]))
    print_stats()
    
    # gamma plot
    plt.figure()
    plt.plot(g_bin_centers, g_tof_hist)
    plt.plot(g_bin_centers, g_hist_fit) 
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
#    plt.savefig(save_dir+dist+'cm_g_fit.png',dpi=500)

    # neutron plot    
    plt.figure()
    plt.plot(n_bin_centers, n_tof_hist)
    plt.plot(n_bin_centers, n_hist_fit, label='convo')
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
#    plt.savefig(save_dir+dist+'cm_n_git.png',dpi=500)

    # full
    plt.figure()
    # scale
    n_full_hist_fit = [n*max(tof_hist)/max(n_full_hist_fit) for n in n_full_hist_fit]
    g_full_hist_fit = [g*max(tof_hist[len(tof_hist)*3./5.:])/max(g_full_hist_fit) for g in g_full_hist_fit]
    plt.plot(bin_centers,tof_hist)
    plt.plot(bin_centers, n_full_hist_fit, linewidth=2, linestyle='--')
    plt.plot(bin_centers, g_full_hist_fit, linewidth=2, linestyle='--')
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
#    plt.savefig(save_dir+dist+'cm_tof_fits.png',dpi=500)

plt.show()
#pickle.dump( means_stds, open( "peak_fit_params_11mev.p", "wb" ) )
#print '\nparams saved to peak_fit_params_11mev.p'