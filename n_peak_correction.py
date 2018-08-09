''' Code for fitting tof spectra 
    Neutron spectrum - set gauss==True for guassian fit, gauss==False for guassian convoluted with a leading exponential
    Gamma spectrum - fit with triple guassian and linear fit for background
'''


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc,erf
from scipy.optimize import curve_fit
import scipy.optimize as opt
import lmfit
from lmfit import Model
import pickle

def gaus_exp_convo(x, a, mu, sigma, gamma, m, b):
    # original
    #erfc_arg = (sigma/tau - (x-mu)/sigma)/np.sqrt(2.)
    #func = a*sigma/tau * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma/tau)**2 - (x-mu)/tau)) * (1-erf(erfc_arg))
    # from lmfit model
    #gss = gamma*sigma*sigma
    #arg1 = gamma*(center + gss/2.0 - x)
    #arg2 = (center + gss - x)/(s2*sigma)
    #return amplitude*(gamma/2) * exp(arg1) * erfc(arg2)

    func = a*(gamma/2.) * np.exp(-gamma*(mu-x+0.5*gamma*sigma*sigma)) * erfc((mu+gamma*sigma*sigma-x)/(sigma*np.sqrt(2.)))
    line = m*x + b
    return func + line

def multi_gaus_exp(x, *p):
    n = (len(p)-2)/4
    a = p[:n]
    mu = p[n:2*n]
    sigma = p[2*n:3*n]
    gamma = p[3*n:4*n]
    res = sum([a[i]*(gamma[i]/2.) * np.exp(-gamma[i]*(mu[i]-x+0.5*gamma[i]*sigma[i]*sigma[i])) * erfc((mu[i]+gamma[i]*sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.))) for i in xrange(n)])
    line = m*x + b
    return res + line

def gaussian( x, a, mu, sigma, m, b):
    #a, mu, sigma = p
    res =   a * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    line = m*x + b
    return res + line

def multi_gauss_line(x, *p):
    n = (len(p)-2)/3
    a = p[:n]
    mu = p[n:2*n]
    sigma = p[2*n:3*n]
    guass = sum([a[i] * np.exp( - (x - mu[i])**2.0 / (2.0 * sigma[i]**2.0) ) for i in xrange(n)])
    m = p[-2]
    b = p[-1]
    line = m*x + b
    return guass + line

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

def minimize_function(p, x, y):
    model = multi_gaus_exp(x,*p)
    diff = model - y
    return np.dot(diff, diff)  # return SSE

e_dir = 'C:/Users/raweldon/Research/TUNL/git_programs/tof_codes/plots/'
gauss = False # true if guass fit, flase if guass-exp convolution fit
plt_save = False # if true save plots
save_params = False # if true save params to pickle

dists = ['179','276','369']


 # get vaules by inspection with plot_tof_hist.py
if gauss == True:
    fit_type='gauss'
    n_ranges=[[250.,380.],[220,380.],[180,380.]]
else:
    fit_type = 'gauss_exp_conv'
    n_ranges=[[250.,380.],[220,380.],[180,380.]]

amp = 1.0
sigma = 1.0
m = 0.1
b = 70
low_bound_2 = 240
high_bound_2 = 275  
low_bound_3 = 200
high_bound_3 = 250       
n_p0s = [ [1.0,1.0,1.0,1.0,1.0,1.0,287.0,291,295,299,301,303,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.],
          [1.0,1.0,1.0,1.0,1.0,1.0,1.0,250,255,259,263,267,269,271,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.], 
          [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,218,222,226,230,234,238,241,242,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.] ]

n_bounds = [ ((0,0,0,0,0,0,                        280,280,280,280,280,280,0,0,0,0,0,0,-10,-10),
              (50000,50000,50000,50000,50000,50000,310,310,310,310,310,310,2,2,2,2,2,2,10,1000)),
             ((0,0,0,0,0,0,0,low_bound_2,low_bound_2,low_bound_2,low_bound_2,low_bound_2,low_bound_2,low_bound_2,0,0,0,0,0,0,0,-10,-10),
              (50000,50000,50000,50000,50000,50000,5000,high_bound_2,high_bound_2,high_bound_2,high_bound_2,high_bound_2,high_bound_2,high_bound_2,2,2,2,2,2,2,2,10,1000)),
             ((0,0,0,0,0,0,0,0,low_bound_3,low_bound_3,low_bound_3,low_bound_3,low_bound_3,low_bound_3,low_bound_3,low_bound_3,0,0,0,0,0,0,0,0,-10,-10),
              (50000,50000,50000,50000,50000,50000,5000,5000,high_bound_3,high_bound_3,high_bound_3,high_bound_3,high_bound_3,high_bound_3,high_bound_3,high_bound_3,2,2,2,2,2,2,2,2,10,1000)) ]

n_min_p0s = [ [1.0,1.0,1.0,1.0,1.0,1.0,287.0,291,295,299,301,303,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.],
          [1.0,1.0,1.0,1.0,1.0,1.0,1.0,250,255,259,263,267,269,271,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.], 
          [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,218,222,226,230,234,238,241,242,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, 0.1,70.] ]

n_bounds_min = [ ((0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(280,310),(280,310),(280,310),(280,310),(280,310),(280,310),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),
                  (0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(-10,10),(-10,1000)),
                 ((0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(low_bound_2,high_bound_2),(low_bound_2,high_bound_2),(low_bound_2,high_bound_2),(low_bound_2,high_bound_2),
                  (low_bound_2,high_bound_2),(low_bound_2,high_bound_2),(low_bound_2,high_bound_2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),
                  (0,5000),(0,5000),(-10,10),(-10,1000)),
                 ((0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(low_bound_3,high_bound_3),(low_bound_3,high_bound_3),(low_bound_3,high_bound_3),(low_bound_3,high_bound_3),
                   (low_bound_3,high_bound_3),(low_bound_3,high_bound_3),(low_bound_3,high_bound_3),(low_bound_3,high_bound_3),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,2),(0,5000),(0,5000),
                   (0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(0,5000),(-10,10),(-10,1000)) ]

for i in xrange(len(n_p0s)):
    print len(n_min_p0s[i]),len(n_bounds_min[i])

means_stds=[]
for index,dist in enumerate(dists):
    print '\n---------------------------------------------------'
    print '---------------------- '+str(dist)+'cm ----------------------'
    print '---------------------------------------------------'
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = get_range(tof, 200, 450)
    n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
    
    # build hists  
    tof_hist, bin_centers = build_hist(tof, 2e3)
    n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
    
    # neutrons
    p0 = n_p0s[index]
    bounds = n_bounds[index]
    if gauss == True:                                  
        coeff, var_matrix = curve_fit(multi_gauss_line, n_bin_centers, n_tof_hist, p0=p0, bounds=bounds, max_nfev=10000)
        n_hist_fit = multi_gauss_line(n_bin_centers, *coeff)
        n_full_hist_fit = multi_gauss_line(bin_centers, *coeff)
        for c in coeff:
            print c

    else:
        res = opt.minimize(minimize_function, x0=n_min_p0s[index], args=(n_bin_centers,n_tof_hist), method='L-BFGS-B', bounds=n_bounds_min[index])
        print res.x
        coeff = res.x
        n_hist_fit = multi_gaus_exp(n_bin_centers, *coeff)
        n_full_hist_fit = multi_gaus_exp(bin_centers, *coeff)
    
    # full
    plt.figure()
    # scale
    n_full_hist_fit = [n*max(tof_hist)/max(n_full_hist_fit) for n in n_full_hist_fit]
    plt.plot(bin_centers,tof_hist)
    plt.plot(bin_centers, n_full_hist_fit, linewidth=2, linestyle='--')
    l = (len(coeff)-2)/4
    for c in xrange(l):
        plt.plot(bin_centers,gaus_exp_convo(bin_centers,coeff[c],coeff[c+l],coeff[c+2*l],coeff[c+3*l],coeff[-2],coeff[-1]),linewidth=2)

    plt.ylabel('counts')
    plt.xlabel('time (ns)')
    plt.xlim(200,320)
    if plt_save == True:
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_tof_fits.png',dpi=500)

plt.show()
if save_params == True:
    pickle.dump( means_stds, open( "peak_fit_params_11mev_"+fit_type+".p", "wb" ) )
    print '\nparams saved to peak_fit_params_11mev_'+fit_type+'.p'