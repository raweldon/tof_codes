''' Code for fitting tof spectra 
    Neutron spectrum - set gauss==True for guassian fit, gauss==False for guassian convoluted with a leading exponential
    Gamma spectrum - fit with double guassian
'''


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc,erf
from scipy.optimize import curve_fit
import lmfit
from lmfit import Model
import pickle

def gaus_exp_convo(x, a, mu, sigma, tau):
    erfc_arg = (sigma/tau - (x-mu)/sigma)/np.sqrt(2.)
    func = a*sigma/tau * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma/tau)**2 - (x-mu)/tau)) * (1-erf(erfc_arg))
    return func

def triple_gaus_exp_convo(x, a_1, mu_1, sigma_1, tau_1, a_2, mu_2, sigma_2, tau_2, a_3, mu_3, sigma_3, tau_3):
    erfc_arg_1 = (sigma_1/tau_1 - (x-mu_1)/sigma_1)/np.sqrt(2.)
    erfc_arg_2 = (sigma_2/tau_2 - (x-mu_2)/sigma_2)/np.sqrt(2.)
    erfc_arg_3 = (sigma_3/tau_3 - (x-mu_3)/sigma_3)/np.sqrt(2.)
    func = ( a_1*sigma_1/tau_1 * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma_1/tau_1)**2 - (x-mu_1)/tau_1)) * (1-erf(erfc_arg_1)) +
             a_2*sigma_2/tau_2 * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma_2/tau_2)**2 - (x-mu_2)/tau_2)) * (1-erf(erfc_arg_2)) +
             a_3*sigma_3/tau_3 * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma_3/tau_3)**2 - (x-mu_3)/tau_3)) * (1-erf(erfc_arg_3)) )
    return func

def gaussian( x, a, mu, sigma):
    #a, mu, sigma = p
    res =   a * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    return res

def double_gauss(x, *p):
    a1, mu1, sigma1, a2, mu2, sigma2 = p
    res = a1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + a2*np.exp(-(x-mu2)**2/(2.*sigma2**2))
    return res

def triple_gauss_line(x, a1, mu1, sigma1, a2, mu2, sigma2, a3, mu3, sigma3, m, b ):
    trip = a1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + a2*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + a3*np.exp(-(x-mu3)**2/(2.*sigma3**2)) 
    line = m*x + b
    return trip + line

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

save_dir = 'C:/Users/raweldon/Research/TUNL/git_programs/tof_codes/plots/'
gauss = True # true if guass fit, flase if guass-exp convolution fit
plt_save = False # if true save plots
save_params = False # if true save params to pickle

dists = ['179','276','369']

# get vaules by inspection with plot_tof_hist.py
if gauss == True:
    fit_type='gauss'
    n_ranges=[[302.,307.],[269,275.],[240.5,246.]]
    g_ranges=[[360.,430.],[360.,430.],[360.,430]]
else:
    fit_type = 'gauss_exp_conv'
    n_ranges=[[302.,307.],[269,275.],[240.5,246.]]
    g_ranges=[[360.,410.],[360.,404.],[360.,402]]

n_p0s = [[1.0, 302.5, 1.0],[1.0, 270., 1.0],[1.0,241.,1.0]]
g_p0s = [[1.0, 398.5, 0.1, 1.0, 1.0, 401.5, 0.1, 1.0, 1.0, 406.3, 0.1, 1.0],
         [1.0, 395.0, 0.1, 1.0, 1.0, 398.5, 0.1, 1.0, 1.0, 403.3, 0.1, 1.0],
         [1.0, 392.0, 0.1, 1.0, 1.0, 395.3, 0.1, 1.0, 1.0, 400.5, 0.1, 1.0]]    

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
    n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
    g_tof_hist, g_bin_centers = build_hist(g_tof, 1000)
    
    # neutrons
    p0 = n_p0s[index]
    if gauss == True:
        gmodel = lmfit.Model(gaussian)
        params = gmodel.make_params(a=1000,mu=p0[1],sigma=p0[2])
    
        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting neutron peak with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=True)
    
        coeff = (res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)
        n_hist_fit = gaussian(n_bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)
        n_full_hist_fit = gaussian(bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)

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

    means_stds.append((coeff[1],coeff[2]))

    # gammas
    p0 = g_p0s[index]
 
    gmodel = Model(triple_gauss_line)
    gmodel.set_param_hint('mu1',value=p0[1], min=392., max=401)
    gmodel.set_param_hint('mu2', value=p0[5], min=400, max=403)
    gmodel.set_param_hint('mu3', value=p0[9], min=405, max=409)
    params = gmodel.make_params(a1=p0[0], mu1=p0[1], sigma1=p0[2], a2=p0[4], mu2=p0[5], sigma2=p0[6],a3=p0[8], mu3=p0[9], sigma3=p0[10], m=1.0, b=1.0)

    res = gmodel.fit(g_tof_hist, params, x=g_bin_centers, nan_policy='omit')
    print '\nFitting gamma peaks with LMFIT'
    print res.message
    print lmfit.fit_report(res,show_correl=True)

    coeff = ( res.params['a1'].value, res.params['mu1'].value, res.params['sigma1'].value, res.params['a2'].value, res.params['mu2'].value, res.params['sigma2'].value,
              res.params['a3'].value, res.params['mu3'].value, res.params['sigma3'].value, res.params['m'].value, res.params['b'].value )
    g_hist_fit = triple_gauss_line( g_bin_centers, res.params['a1'].value, res.params['mu1'].value, res.params['sigma1'].value, 
                                                       res.params['a2'].value, res.params['mu2'].value, res.params['sigma2'].value,
                                                       res.params['a3'].value, res.params['mu3'].value, res.params['sigma3'].value, res.params['m'].value, res.params['b'].value )
    g_full_hist_fit = triple_gauss_line( bin_centers, res.params['a1'].value, res.params['mu1'].value, res.params['sigma1'].value, 
                                                          res.params['a2'].value, res.params['mu2'].value, res.params['sigma2'].value,
                                                          res.params['a3'].value, res.params['mu3'].value, res.params['sigma3'].value, res.params['m'].value, res.params['b'].value )
 
#    gmodel = lmfit.Model(triple_gaus_exp_convo)
#    params = gmodel.make_params(a_1=p0[0], mu_1=p0[1], sigma_1=p0[2], tau_1=p0[3], a_2=p0[4], mu_2=p0[5], sigma_2=p0[6], tau_2=p0[7],
#                                a_3=p0[8], mu_3=p0[9], sigma_3=p0[10], tau_3=p0[11])
#
#    res = gmodel.fit(g_tof_hist, params, x=g_bin_centers, nan_policy='omit')
#    print '\nFitting gamma peaks with LMFIT'
#    print res.message
#    print lmfit.fit_report(res,show_correl=True)
#
#    coeff = ( res.params['a_1'].value, res.params['mu_1'].value, res.params['sigma_1'].value, res.params['tau_1'].value,
#              res.params['a_2'].value, res.params['mu_2'].value, res.params['sigma_2'].value, res.params['tau_2'].value,
#              res.params['a_3'].value, res.params['mu_3'].value, res.params['sigma_3'].value, res.params['tau_3'].value )
#    g_hist_fit = triple_gaus_exp_convo( g_bin_centers, res.params['a_1'].value, res.params['mu_1'].value, res.params['sigma_1'].value, res.params['tau_1'].value,
#                                        res.params['a_2'].value, res.params['mu_2'].value, res.params['sigma_2'].value, res.params['tau_2'].value,
#                                        res.params['a_3'].value, res.params['mu_3'].value, res.params['sigma_3'].value, res.params['tau_3'].value  )
#    g_full_hist_fit = triple_gaus_exp_convo( bin_centers, res.params['a_1'].value, res.params['mu_1'].value, res.params['sigma_1'].value, res.params['tau_1'].value,
#                                             res.params['a_2'].value, res.params['mu_2'].value, res.params['sigma_2'].value, res.params['tau_2'].value,
#                                             res.params['a_3'].value, res.params['mu_3'].value, res.params['sigma_3'].value, res.params['tau_3'].value )
    means_stds.append((coeff[1],coeff[2],coeff[4],coeff[5],coeff[7],coeff[8]))
    
    # gamma plot
    fig1, ax1 = plt.subplots()
    plt.plot(g_bin_centers, g_tof_hist)
    plt.plot(g_bin_centers, g_hist_fit) 
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
    plt.text(0.7,0.75,'$\mu_1 = $'+str(round(means_stds[1+2*index][0],3))+'\n$\sigma_1 =$ '+str(round(means_stds[1+2*index][1],3))+'\n$\mu_2 =$ '+
             str(round(means_stds[1+2*index][2],3))+'\n$\sigma_2 =$ '+str(round(means_stds[1+2*index][3],3))+'\n$\mu_3 =$ '+
             str(round(means_stds[1+2*index][4],3))+'\n$\sigma_3 =$ '+str(round(means_stds[1+2*index][5],3)),transform=ax1.transAxes)
    if plt_save == True:
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_g_fit.png',dpi=500)

    # neutron plot    
    fig2, ax2 = plt.subplots()
    plt.plot(n_bin_centers, n_tof_hist)
    plt.plot(n_bin_centers, n_hist_fit, label='convo')
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
    plt.text(0.7,0.75,'$\mu =$ '+str(round(means_stds[0+2*index][0],3))+'\n$\sigma =$ '+str(round(means_stds[0+2*index][1],3)), transform=ax2.transAxes)
    if plt_save == True:    
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_n_fit.png',dpi=500)

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
    if plt_save == True:
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_tof_fits.png',dpi=500)

plt.show()
if save_params == True:
    pickle.dump( means_stds, open( "peak_fit_params_11mev_"+fit_type+".p", "wb" ) )
    print '\nparams saved to peak_fit_params_11mev_'+fit_type+'.p'