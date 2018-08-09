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
import collections

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

def gaussian( x, a, mu, sigma, m, b):
    #a, mu, sigma = p
    res =   a * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    line = m*x + b
    return res + line

def gaus_6(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3, a4, mu4, sigma4, a5, mu5, sigma5, a6, mu6, sigma6, m, b ):
    a = [a1,a2,a3,a4,a5,a6]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6]
    func = sum([a[i] * np.exp( - (x - mu[i])**2.0 / (2.0 * sigma[i]**2.0) ) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_7(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3, a4, mu4, sigma4, a5, mu5, sigma5, a6, mu6, sigma6, a7, mu7, sigma7, m, b ):
    a = [a1,a2,a3,a4,a5,a6,a7]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6,mu7]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7]
    func = sum([a[i] * np.exp( - (x - mu[i])**2.0 / (2.0 * sigma[i]**2.0) ) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_8(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3, a4, mu4, sigma4, a5, mu5, sigma5, a6, mu6, sigma6, 
           a7, mu7, sigma7, a8, mu8, sigma8, m, b ):
    a = [a1,a2,a3,a4,a5,a6,a7,a8]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8]
    func = sum([a[i] * np.exp( - (x - mu[i])**2.0 / (2.0 * sigma[i]**2.0) ) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_exp_6(x, a1, mu1, sigma1, gamma1, a2, mu2, sigma2, gamma2, a3, mu3, sigma3, gamma3, 
               a4, mu4, sigma4, gamma4, a5, mu5, sigma5, gamma5, a6, mu6, sigma6, gamma6, m, b ):
    a = [a1,a2,a3,a4,a5,a6]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6]
    gamma = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6]
    func = [a[i]*(gamma[i]/2.) * np.exp(-gamma[i]*(mu[i]-x+0.5*gamma[i]*sigma[i]*sigma[i])) * 
            erfc((mu[i]+gamma[i]*sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.))) for i in xrange(len(a))]
    line = m*x + b
    return func + line

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

e_dir = 'C:/Users/raweldon/Research/TUNL/git_programs/tof_codes/plots/'
gauss = True # true if guass fit, flase if guass-exp convolution fit
plt_save = False # if true save plots
save_params = False # if true save params to pickle

dists = ['179','276','369']

 # get vaules by inspection with plot_tof_hist.py
if gauss == True:
    fit_type='gauss'
    n_ranges=[[255.,380.],[205,380.],[160,380.]]
else:
    fit_type = 'gauss_exp_conv'
    n_ranges=[[250.,380.],[220,380.],[180,380.]]
     
n_p0s = [ [1.0, 287.0, 1.0,
           1.0, 291.0, 1.0,
           1.0, 295.0, 1.0,
           1.0, 299.0, 1.0,
           1.0, 301.0, 1.0,
           1.0 ,303.0 ,1.0, 0.1,100.],
          [1.0, 251.0, 1.0,
           1.0, 255.0, 1.0,
           1.0, 259.0, 1.0,
           1.0, 263.0, 1.0,
           1.0, 267.0, 1.0,
           1.0, 269.0, 1.0,
           1.0, 271.0, 1.0, 0.1,100.], 
          [1.0, 218.0, 1.0,
           1.0, 222.0, 1.0,
           1.0, 226.0, 1.0,
           1.0, 230.0, 1.0,
           1.0, 234.0, 1.0,
           1.0, 238.0, 1.0,
           1.0, 241.0, 1.0,
           1.0, 242.0, 1.0, 0.1,100] ]

n_bounds = [ [(0,50000), (280,310), (0,2),
              (0,50000), (280,310), (0,2),
              (0,50000), (280,310), (0,2),
              (0,50000), (280,310), (0,2),
              (0,50000), (280,310), (0,2),
              (0,50000), (280,310), (0,2), (-0.04,0.0),(-10,300)],
             [(0,50000), (240,275), (0,2),
              (0,50000), (240,275), (0,2),
              (0,50000), (240,275), (0,2),
              (0,50000), (240,275), (0,2),
              (0,50000), (240,275), (0,2),
              (0,50000), (240,275), (0,2), 
              (0,50000), (240,275), (0,2), (-0.04,0.0),(-10,200)],
             [(0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2), 
              (0,50000), (200,250), (0,2),
              (0,50000), (200,250), (0,2), (-0.04,0.0),(-10,300)] ] 

means_stds=[]
model = (gaus_6,gaus_7,gaus_8)
for index,dist in enumerate(dists):
    print '\n---------------------------------------------------'
    print '---------------------- '+str(dist)+'cm ----------------------'
    print '---------------------------------------------------'
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = get_range(tof, 150, 450)
    n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
    
    # build hists  
    tof_hist, bin_centers = build_hist(tof, 2e3)
    n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
    
    # neutrons
    p0 = n_p0s[index]
    bounds = n_bounds[index]
    if gauss == True:
        gmodel = lmfit.Model(model[index])
        gmodel.set_param_hint('a1', value=p0[0], min=bounds[0][0], max=bounds[0][1])
        gmodel.set_param_hint('mu1', value=p0[1], min=bounds[1][0], max=bounds[1][1])
        gmodel.set_param_hint('sigma1', value=p0[2], min=bounds[2][0], max=bounds[2][1])
        gmodel.set_param_hint('a2', value=p0[3], min=bounds[3][0], max=bounds[3][1])
        gmodel.set_param_hint('mu2', value=p0[4], min=bounds[4][0], max=bounds[4][1])
        gmodel.set_param_hint('sigma2', value=p0[5], min=bounds[5][0], max=bounds[5][1])
        gmodel.set_param_hint('a3', value=p0[6], min=bounds[6][0], max=bounds[6][1])
        gmodel.set_param_hint('mu3', value=p0[7], min=bounds[7][0], max=bounds[7][1])
        gmodel.set_param_hint('sigma3', value=p0[8], min=bounds[8][0], max=bounds[8][1])
        gmodel.set_param_hint('a4', value=p0[9], min=bounds[9][0], max=bounds[9][1])
        gmodel.set_param_hint('mu4', value=p0[10], min=bounds[10][0], max=bounds[10][1])
        gmodel.set_param_hint('sigma4', value=p0[11], min=bounds[11][0], max=bounds[11][1])
        gmodel.set_param_hint('a5', value=p0[12], min=bounds[12][0], max=bounds[12][1])
        gmodel.set_param_hint('mu5', value=p0[13], min=bounds[13][0], max=bounds[13][1])
        gmodel.set_param_hint('sigma5', value=p0[14], min=bounds[14][0], max=bounds[14][1])
        gmodel.set_param_hint('a6', value=p0[15], min=bounds[15][0], max=bounds[15][1])
        gmodel.set_param_hint('mu6', value=p0[16], min=bounds[16][0], max=bounds[16][1])
        gmodel.set_param_hint('sigma6', value=p0[17], min=bounds[17][0], max=bounds[17][1])

        if index == 0:
            gmodel.set_param_hint('m', value=p0[18], min=bounds[18][0], max=bounds[18][1])
            gmodel.set_param_hint('b', value=p0[19], min=bounds[19][0], max=bounds[19][1])
        if index == 1:
            gmodel.set_param_hint('a7', value=p0[18], min=bounds[18][0], max=bounds[18][1])
            gmodel.set_param_hint('mu7', value=p0[19], min=bounds[19][0], max=bounds[19][1])
            gmodel.set_param_hint('sigma7', value=p0[20], min=bounds[20][0], max=bounds[20][1])
            gmodel.set_param_hint('m', value=p0[21], min=bounds[21][0], max=bounds[21][1])
            gmodel.set_param_hint('b', value=p0[22], min=bounds[22][0], max=bounds[22][1])
        if index == 2:
            gmodel.set_param_hint('a7', value=p0[18], min=bounds[18][0], max=bounds[18][1])
            gmodel.set_param_hint('mu7', value=p0[19], min=bounds[19][0], max=bounds[19][1])
            gmodel.set_param_hint('sigma7', value=p0[20], min=bounds[20][0], max=bounds[20][1])
            gmodel.set_param_hint('a8', value=p0[21], min=bounds[21][0], max=bounds[21][1])
            gmodel.set_param_hint('mu8', value=p0[22], min=bounds[22][0], max=bounds[22][1])
            gmodel.set_param_hint('sigma8', value=p0[23], min=bounds[23][0], max=bounds[23][1])
            gmodel.set_param_hint('m', value=p0[24], min=bounds[24][0], max=bounds[24][1])
            gmodel.set_param_hint('b', value=p0[25], min=bounds[25][0], max=bounds[25][1])

        params = gmodel.make_params()
    
        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting neutron peak with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=False)
        
        f = res.best_values

        plt.figure()   
        plt.plot(n_bin_centers,n_tof_hist)
        plt.plot(n_bin_centers,res.best_fit)
        for i in xrange(1,7):
            plt.plot(n_bin_centers,gaussian(n_bin_centers,f['a'+str(i)],f['mu'+str(i)],f['sigma'+str(i)],f['m'],f['b']))
    

    else:
        gmodel = lmfit.Model(gaus_exp_convo)
        params = gmodel.make_params(a=1000,mu=p0[1],sigma=p0[2],tau=p0[3])

        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting neutron peak with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=True)

        coeff = (res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
        n_hist_fit = gaus_exp_convo(n_bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
        n_full_hist_fit = gaus_exp_convo(bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
plt.show()

    # full
#    plt.figure()
#    # scale
#    n_full_hist_fit = [n*max(tof_hist)/max(n_full_hist_fit) for n in n_full_hist_fit]
#    plt.plot(bin_centers,tof_hist)
#    plt.plot(bin_centers, n_full_hist_fit, linewidth=2, linestyle='--')
#    l = (len(coeff)-2)/4
#    for c in xrange(l):
#        plt.plot(bin_centers,gaus_exp_convo(bin_centers,coeff[c],coeff[c+l],coeff[c+2*l],coeff[c+3*l],coeff[-2],coeff[-1]),linewidth=2)
#
#    plt.ylabel('counts')
#    plt.xlabel('time (ns)')
#    plt.xlim(200,320)
#    if plt_save == True:
#        plt.savefig(save_dir+fit_type+'_'+dist+'cm_tof_fits.png',dpi=500)
#
#plt.show()
#if save_params == True:
#    pickle.dump( means_stds, open( "peak_fit_params_11mev_"+fit_type+".p", "wb" ) )
#    print '\nparams saved to peak_fit_params_11mev_'+fit_type+'.p'