''' New code for tof appendix of dissertation
    multiple gaussia exponential convolution fits
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

    func = a*(gamma/2.) * np.exp(-gamma*(mu-x+0.5*gamma*sigma*sigma)) * (1 - erf(-(mu+gamma*sigma*sigma-x)/(sigma*np.sqrt(2.))))
    line = m*x + b
    return func  + line


def gaus_exp_2(x, a1, mu1, sigma1, gamma, a2, mu2, sigma2, m, b ):
    a = [a1,a2]
    mu = [mu1,mu2]
    sigma = [sigma1,sigma2]
    func = sum([a[i]*(gamma/2.) * np.exp(-gamma*(mu[i]-x+0.5*gamma*sigma[i]*sigma[i])) * 
            (1 - erf(-(mu[i]+sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.)))) for i in xrange(len(a))])
    line = m*x + b
    return func# + line

def gaus_exp_5(x, a1, mu1, sigma1, gamma, a2, mu2, sigma2, a3, mu3, sigma3, 
               a4, mu4, sigma4, a5, mu5, sigma5, m, b ):
    a = [a1,a2,a3,a4,a5]
    mu = [mu1,mu2,mu3,mu4,mu5]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5]
    func = sum([a[i]*(gamma/2.) * np.exp(-gamma*(mu[i]-x+0.5*gamma*sigma[i]*sigma[i])) * 
            (1 - erf(-(mu[i]+sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.)))) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_exp_6(x, a1, mu1, sigma1, gamma, a2, mu2, sigma2, a3, mu3, sigma3, 
               a4, mu4, sigma4, a5, mu5, sigma5, a6, mu6, sigma6, m, b ):
    a = [a1,a2,a3,a4,a5,a6]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6]
    func = sum([a[i]*(gamma/2.) * np.exp(-gamma*(mu[i]-x+0.5*gamma*sigma[i]*sigma[i])) * 
            (1 - erf(-(mu[i]+sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.)))) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_exp_7(x, a1, mu1, sigma1, gamma, a2, mu2, sigma2, a3, mu3, sigma3,  
               a4, mu4, sigma4, a5, mu5, sigma5, a6, mu6, sigma6, a7, mu7, sigma7, m, b ):
    a = [a1,a2,a3,a4,a5,a6,a7]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6,mu7]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7]
    func = sum([a[i]*(gamma/2.) * np.exp(-gamma*(mu[i]-x+0.5*gamma*sigma[i]*sigma[i])) * 
            erfc((mu[i]+gamma*sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.))) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def gaus_exp_8(x, a1, mu1, sigma1, gamma, a2, mu2, sigma2, a3, mu3, sigma3, 
               a4, mu4, sigma4, a5, mu5, sigma5,a6, mu6, sigma6, 
               a7, mu7, sigma7, a8, mu8, sigma8, m, b ):
    a = [a1,a2,a3,a4,a5,a6,a7,a8]
    mu = [mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8]
    sigma = [sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8]
    func = sum([a[i]*(gamma/2.) * np.exp(-gamma*(mu[i]-x+0.5*gamma*sigma[i]*sigma[i])) * 
            erfc((mu[i]+gamma*sigma[i]*sigma[i]-x)/(sigma[i]*np.sqrt(2.))) for i in xrange(len(a))])
    line = m*x + b
    return func + line

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

def multi_fit(n_p0s, n_bounds, n_ranges, peaks, model, dists, e_dir, plt_save, save_params):

    max_vals = []
    for index, dist in enumerate(dists):
        print '\n---------------------------------------------------'
        print '---------------------- '+str(dist)+'cm ----------------------'
        print '---------------------------------------------------'
        tof_spec = np.load('dist_'+dist+'.npz')
        tof = tof_spec['data']
        tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
        tof = get_range(tof, 150, 450)
        n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
        
        # build hists  
        tof_hist, bin_centers = build_hist(tof, 2000)
        n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
        
        # neutrons
        p0 = n_p0s[index]
        bounds = n_bounds[index]
    
        gmodel = lmfit.Model(model[index])
        gmodel.set_param_hint('a1', value=p0[0], min=bounds[0][0], max=bounds[0][1])
        gmodel.set_param_hint('mu1', value=p0[1], min=bounds[1][0], max=bounds[1][1])
        gmodel.set_param_hint('sigma1', value=p0[2], min=bounds[2][0], max=bounds[2][1],    vary=True)
        gmodel.set_param_hint('gamma', value=p0[3], min=bounds[3][0], max=bounds[3][1],    vary=True)
        gmodel.set_param_hint('a2', value=p0[4], min=bounds[4][0], max=bounds[4][1])
        gmodel.set_param_hint('mu2', value=p0[5], min=bounds[5][0], max=bounds[5][1])
        gmodel.set_param_hint('sigma2', value=p0[6], min=bounds[6][0], max=bounds[6][1],    vary=True)
        gmodel.set_param_hint('a3', value=p0[8], min=bounds[8][0], max=bounds[8][1])
        gmodel.set_param_hint('mu3', value=p0[9], min=bounds[9][0], max=bounds[9][1])
        gmodel.set_param_hint('sigma3', value=p0[10], min=bounds[10][0], max=bounds[10][1], vary=True)
        gmodel.set_param_hint('a4', value=p0[12], min=bounds[12][0], max=bounds[12][1])
        gmodel.set_param_hint('mu4', value=p0[13], min=bounds[13][0], max=bounds[13][1])
        gmodel.set_param_hint('sigma4', value=p0[14], min=bounds[14][0], max=bounds[14][1], vary=True)
        gmodel.set_param_hint('a5', value=p0[16], min=bounds[16][0], max=bounds[16][1])
        gmodel.set_param_hint('mu5', value=p0[17], min=bounds[17][0], max=bounds[17][1])
        gmodel.set_param_hint('sigma5', value=p0[18], min=bounds[18][0], max=bounds[18][1], vary=True)

        gmodel.set_param_hint('a6', value=p0[20], min=bounds[20][0], max=bounds[20][1])
        gmodel.set_param_hint('mu6', value=p0[21], min=bounds[21][0], max=bounds[21][1])
        gmodel.set_param_hint('sigma6', value=p0[22], min=bounds[22][0], max=bounds[22][1], vary=True)
        gmodel.set_param_hint('m', value=p0[24], min=bounds[24][0], max=bounds[24][1])
        gmodel.set_param_hint('b', value=p0[25], min=bounds[25][0], max=bounds[25][1])

        params = gmodel.make_params()

        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting neutron peak with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=False)
        f = res.best_values

        plt.figure()   
        plt.plot(n_bin_centers, n_tof_hist, alpha=0.9, label='measured spectrum')
        plt.plot(n_bin_centers, res.best_fit, '--', linewidth=1.5, zorder=300, label='fit')
        for i in xrange(1,peaks[index]+1):
            plt.plot(n_bin_centers,gaus_exp_convo(n_bin_centers,f['a'+str(i)],f['mu'+str(i)],f['sigma'+str(i)],f['gamma'],f['m'],f['b']), alpha=0.9)
        #plt.xlim(270, 315)
        plt.ylim(0, 1.3e4)
        plt.xlabel('Time (ns)', fontsize=18)
        plt.ylabel('Counts', fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        plt.tight_layout()
        #plt.savefig('n_peak_fit_179cm.pdf')

        # print max values of convolution
        xvals = np.linspace(230, 310, 1000)
        max_idx = np.argmax(gaus_exp_convo(xvals , f['a6'], f['mu6'], f['sigma6'], f['gamma'], f['m'], f['b']))
        max_vals.append(xvals[max_idx])

    if len(dists) == 3:
        print '\n  179 cm       276 cm       369 cm'
        print '{:^8.2f} ns {:>8.2f} ns {:>8.2f} ns'. format(max_vals[0], max_vals[1], max_vals[2])
    plt.show()


def zoomed_fits(n_p0s, n_bounds, n_ranges, peaks, model, dists, e_dir, plt_save, save_params):
    for index, dist in enumerate(dists):
        print '\n---------------------------------------------------'
        print '---------------------- '+str(dist)+'cm ----------------------'
        print '---------------------------------------------------'
        tof_spec = np.load('dist_'+dist+'.npz')
        tof = tof_spec['data']
        tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
        tof = get_range(tof, 150, 450)
        n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
        
        # build hists  
        tof_hist, bin_centers = build_hist(tof, 2000)
        n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
        
        # neutrons
        p0 = n_p0s[index]
        bounds = n_bounds[index]
    
        gmodel = lmfit.Model(model[index])
        gmodel.set_param_hint('a1', value=p0[0], min=bounds[0][0], max=bounds[0][1])
        gmodel.set_param_hint('mu1', value=p0[1], min=bounds[1][0], max=bounds[1][1])
        gmodel.set_param_hint('sigma1', value=p0[2], min=bounds[2][0], max=bounds[2][1],    vary=True)
        gmodel.set_param_hint('gamma', value=p0[3], min=bounds[3][0], max=bounds[3][1],    vary=True)
        gmodel.set_param_hint('a2', value=p0[4], min=bounds[4][0], max=bounds[4][1])
        gmodel.set_param_hint('mu2', value=p0[5], min=bounds[5][0], max=bounds[5][1])
        gmodel.set_param_hint('sigma2', value=p0[6], min=bounds[6][0], max=bounds[6][1],    vary=True)
        #gmodel.set_param_hint('m', value=p0[8], min=bounds[8][0], max=bounds[8][1])
        #gmodel.set_param_hint('b', value=p0[9], min=bounds[9][0], max=bounds[9][1])


        params = gmodel.make_params()

        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting neutron peak with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=False)
        f = res.best_values

        plt.figure()   
        plt.plot(n_bin_centers,n_tof_hist)
        plt.plot(n_bin_centers,res.best_fit,'--',linewidth=2)
        for i in xrange(1,peaks[index]+1):
            plt.plot(n_bin_centers,gaus_exp_convo(n_bin_centers,f['a'+str(i)],f['mu'+str(i)],f['sigma'+str(i)],f['gamma']))#,f['m'],f['b']))
        plt.ylim(0, 14e3)

    plt.show()

if __name__ == '__main__':

    e_dir = 'C:/Users/raweldon/Research/TUNL/git_programs/tof_codes/plots/'
    plt_save = False # if true save plots
    save_params = False # if true save params to pickle

    dists = ['179','276','369']
    #dists = ['179']

    # get vaules by inspection with plot_tof_hist.py
    model = (gaus_exp_6, gaus_exp_6, gaus_exp_6)
    peaks = (6, 6, 6)
    n_ranges=[[260., 320.],[225., 300.],[160., 300.]]

    n_p0s = np.array([ [1000.0, 285.0, 1.1, 0.31, 
                        1000.0, 289.0, 1.0, 0.31,
                        1000.0, 293.0, 1.0, 0.31,
                        1000.0, 297.0, 1.0, 0.31,
                        1000.0, 301.0, 1.0, 0.31,
                        1000.0, 303.0 ,1.0, 0.31, -0.1, 100.],
                       [1000.0, 254.0, 1.0, 0.31,
                        1000.0, 258.0, 1.0, 0.31,
                        1000.0, 262.0, 1.0, 0.31,
                        1000.0, 266.0, 1.0, 0.31,
                        1000.0, 270.0, 1.0, 0.31,
                        1000.0, 272.0, 1.0, 0.31, -0.1, 100.], 
                       [1000.0, 222.0, 1.0, 0.31,
                        1000.0, 226.0, 1.0, 0.31,
                        1000.0, 230.0, 1.0, 0.31,
                        1000.0, 234.0, 1.0, 0.31,
                        1000.0, 238.0, 1.0, 0.31,
                        1000.0, 241.0, 1.0, 0.31, -0.1, 100.] ], dtype=np.object)

    n_bounds = [ [(0, 200000), (280, 290), (0.5, 3), (0.3, 0.5), 
                  (0, 200000), (284, 294), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (288, 298), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (292, 302), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (396, 304), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (302, 305), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 300)],
                 [(0, 200000), (240, 270), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (240, 270), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (240, 270), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (240, 270), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (266, 271), (0.5, 3), (0.1, 0.5), 
                  (0, 200000), (270, 275), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 200)],
                 [(0, 200000), (210, 228), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (222, 230), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (226, 240), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (230, 240), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (238, 241), (0.5, 3), (0.1, 0.5),
                  (0, 200000), (240, 242), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 300)] ]

    multi_fit(n_p0s, n_bounds, n_ranges, peaks, model, dists, e_dir, plt_save, save_params)

    # n_p0s = np.array([ [1000.0, 301.0, 1.0, 0.31,
    #                     1000.0, 303.0 ,1.0, 0.31, -0.1, 100.],
    #                    [1000.0, 270.0, 1.0, 0.31,
    #                     1000.0, 272.0, 1.0, 0.31, -0.1, 100.], 
    #                    [1000.0, 238.0, 1.0, 0.31,
    #                     1000.0, 241.0, 1.0, 0.31, -0.1, 100.] ], dtype=np.object)

    # n_bounds = [ [(0, 200000), (296, 304), (0.5, 3), (0.1, 0.5),
    #               (0, 200000), (302, 305), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 300)],
    #              [(0, 200000), (266, 276), (0.5, 3), (0.1, 0.5),
    #               (0, 200000), (270, 273), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 200)],
    #              [(0, 200000), (234, 244), (0.5, 3), (0.1, 0.5),
    #               (0, 200000), (240, 243), (0.5, 3), (0.1, 0.5), (-0.1, 0.0), (-10, 300)] ]

    # zoomed_model = (gaus_exp_2, gaus_exp_2, gaus_exp_2)
    # zoomed_peaks = (2, 2, 2)
    # zoomed_ranges= ((295, 310), (264, 279), (232, 250))
    # zoomed_fits(n_p0s, n_bounds, n_ranges, zoomed_peaks, zoomed_model, dists, e_dir, plt_save, save_params)
    