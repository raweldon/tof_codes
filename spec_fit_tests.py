import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.optimize import curve_fit

def gaus_exp_convo(x, *p):
    '''not working so well...
       tried reducing fit area of both neutron and gamma peaks, but it tends more toward the exponential
       and does not fit the gaussian at the end well
    '''
    a, tau, mu, sigma = p
    x_erf = (mu - tau*sigma**2 - x)/(np.sqrt(2.)*sigma)
    func = a*np.exp(-tau*(mu - tau*sigma**2/2. - x))/2. * (1. + sp.erf(x_erf))
    return func

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
    print "\nneutron gaussian fit:\n  mu = "+str(coeff[2])+" ns\n  sigma = "+str(coeff[3]
            )+" ns"+'\n  tau = '+str(coeff[1])+'\n amp = '+str(coeff[0])    

#dists = ['175', '235', '284']
dists = ['284']
for dist in dists:
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = get_range(tof, 200, 450)
    n_tof = get_range(tof, 304, 310)
    g_tof = get_range(tof, 384.5, 390)
    
    # build hists  
    tof_hist, bin_centers = build_hist(tof, 5e3)
    n_tof_hist, n_bin_centers = build_hist(n_tof, 300)
    g_tof_hist, g_bin_centers = build_hist(g_tof, 300)
    
    # neutrons
    p0 = [1.0, 305., 1.0] 
    coeff, var_matrix = curve_fit(gaussian, n_bin_centers, n_tof_hist, p0=p0)
    n_hist_fit = gaussian(n_bin_centers, *coeff)
    n_full_hist_fit = gaussian(bin_centers, *coeff)
#    print_stats()

    p0 = [1.0, 1.0, 320., 3.0] 
    coeff, var_matrix = curve_fit(gaus_exp_convo, bin_centers, tof_hist, p0=p0)
    tof_hist_fit = gaus_exp_convo(bin_centers, *coeff)
    print_stats()

    # gamma
    p0 = [1.0, 385., 0.1, 1.0, 388., 0.1] 
    coeff, var_matrix = curve_fit(double_gauss, g_bin_centers, g_tof_hist, p0=p0)
    g_hist_fit = double_gauss(g_bin_centers, *coeff)
    g_full_hist_fit = double_gauss(bin_centers, *coeff)
    print_stats()
    
    plt.figure()
    plt.plot(n_bin_centers, n_tof_hist)
    plt.plot(n_bin_centers, n_hist_fit, label='convo')
    plt.legend()

    plt.figure()
    # scale
    n_full_hist_fit = [n*max(tof_hist)/max(n_full_hist_fit) for n in n_full_hist_fit]
    g_full_hist_fit = [g*max(tof_hist[3000:])/max(g_full_hist_fit) for g in g_full_hist_fit]
    plt.plot(bin_centers,tof_hist)
    plt.plot(bin_centers, n_full_hist_fit, linewidth=2, linestyle='--')
    plt.plot(bin_centers, g_full_hist_fit, linewidth=2, linestyle='--')
    
    # gamma
    plt.figure()
    plt.plot(g_bin_centers, g_tof_hist)
    plt.plot(g_bin_centers, g_hist_fit)    
#    
#    plt.figure()
#    plt.plot(g_short_bin_centers, g_short_hist)
#    plt.plot(g_short_bin_centers, g_short_hist_fit)    


#%% original tries - bad fits
#    # build hists  
#    tof_hist, bin_centers = build_hist(tof, 5e3)
#    n_tof_hist, n_bin_centers = build_hist(n_tof, 300)
#    g1_tof_hist, g1_bin_centers = build_hist(g1_tof, 300)
#    g_short_hist, g_short_bin_centers = build_hist(g_short, 300)    
#    
#    # fit gamma peaks with convolution of gaussian and exponential
#    p0 = [1.0, 0.01, 305., 1.0] 
#    coeff, var_matrix = curve_fit(gaus_exp_convo, n_bin_centers, n_tof_hist, p0=p0)
#    n_hist_fit = gaus_exp_convo(n_bin_centers, *coeff)
#    print_stats()
#
#    p0 = [1.0, 1.0, 320., 3.0] 
#    coeff, var_matrix = curve_fit(gaus_exp_convo, bin_centers, tof_hist, p0=p0)
#    tof_hist_fit = gaus_exp_convo(bin_centers, *coeff)
#    print_stats()
#
#    # gamma
#    p0 = [1.0, 1.0, 385., 1.0] 
#    coeff, var_matrix = curve_fit(gaus_exp_convo, g1_bin_centers, g1_tof_hist, p0=p0)
#    g1_hist_fit = gaus_exp_convo(g1_bin_centers, *coeff)
#    print_stats()
#    
#    p0 = [1.0, 0.1, 385., 1.0] 
#    coeff, var_matrix = curve_fit(gaus_exp_convo, g_short_bin_centers, g_short_hist, p0=p0)
#    g_short_hist_fit = gaus_exp_convo(g_short_bin_centers, *coeff)
#    print_stats()
#     
#    plt.figure()
#    plt.plot(bin_centers,tof_hist)
#    plt.plot(bin_centers, tof_hist_fit)
##    plt.plot(bin_centers,gaus_exp_convo(bin_centers, 10e5, 5.0, 320., 2.0), linewidth=2)
#
#    plt.figure()
#    plt.plot(n_bin_centers, n_tof_hist)
#    plt.plot(n_bin_centers, n_hist_fit, label='convo')
#    plt.legend()
#
#    # gamma
#    plt.figure()
#    plt.plot(g1_bin_centers, g1_tof_hist)
#    plt.plot(g1_bin_centers, g1_hist_fit)    
#    
#    plt.figure()
#    plt.plot(g_short_bin_centers, g_short_hist)
#    plt.plot(g_short_bin_centers, g_short_hist_fit)