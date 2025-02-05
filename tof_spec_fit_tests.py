''' Code for test multi gaussian fit to gamma peaks
    Havar and beam stop peaks are overlapping -- didn't really get it to work, needed to move on
'''


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc,erf
from scipy.optimize import curve_fit
import lmfit
import pickle

def gaus_exp_convo(x, a, mu, sigma, tau):
    erfc_arg = (sigma/tau - (x-mu)/sigma)/np.sqrt(2.)
    func = a*sigma/tau * np.sqrt(np.pi/2.) * np.exp(-(0.5*(sigma/tau)**2 - (x-mu)/tau)) * (1-erf(erfc_arg))
    return func

def gaussian( x, a, mu, sigma):
    #a, mu, sigma = p
    res =   a * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    return res

def double_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2 = p
    res = A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))
    return res

def quad_gauss(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3, A4, mu4, sigma4 = p
    res = A1*np.exp(-(x-mu1)**2/(2.*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2.*sigma2**2)) + A3*np.exp(-(x-mu3)**2/(2.*sigma3**2)) + A4*np.exp(-(x-mu4)**2/(2.*sigma4**2))
    return res

def get_range(vals, low_val, high_val):
    r = [i for i in vals if i>low_val and i<high_val]
    return r

def build_hist(vals, bin_no):
    hist, bins = np.histogram(vals, bin_no)
    bin_centers = (bins[:-1] + bins[1:])/2
    return hist,bin_centers

def print_stats():
    print "\ngamma gaussian fit:\n  mu1 = "+str(coeff[1])+" ns\n  sigma1 = "+str(coeff[2])+" ns"+'\n  mu2 = '+str(coeff[4])+'\n  sigma2 = '+str(coeff[5]
            ) + "  \nmu3 = "+str(coeff[7])+" ns\n  sigma1 = "+str(coeff[8])+" ns"+'\n  mu4 = '+str(coeff[10])+'\n  sigma4 = '+str(coeff[11])

save_dir = 'C:/Users/raweldon/Research/TUNL/git_programs/tof_codes/plots/'
gauss = False # true if guass fit, flase if guass-exp convolution fit
plt_save = False # if true save plots
save_params = False # if true save params to pickle

# 11 MeV
dists = ['180', '256', '363']
# 4 MeV
#dists = ['179','276','369']

# get vaules by inspection with plot_tof_hist.py
# 11.325 MeV
# gauss ranges
if gauss:
    fit_type = 'gauss'
    # values fit top of neutron spectrum
#    n_ranges=[[338.5,341.5],[322.,325.],[298.5,301.5]]
#    g_ranges=[[400.,404.],[397.,402.5],[393.5,398.]]
#    n_p0s = [[1.0, 340., 1.0, 1.0],[1.0, 324., 1.0, 1.0],[1.0, 300., 1.0, 1.0]]
#    g_p0s = [[1.0, 401., 0.1, 1.0, 403.2, 0.1],[1.0, 398.3, 0.1, 1.0, 401., 0.1],[1.0, 394.7, 0.1, 1.0, 397.3, 0.1]]

    # values fit right side of spectrum (early arriving neutrons)
    n_ranges=[[339.7, 345.5],[323.3, 329.],[299.5, 307.5]]
    g_ranges=[[395., 406.],[393., 404],[391, 398.]]
    
    n_p0s = [[1.0, 340., 1.0, 1.0],[1.0, 324., 1.0, 1.0],[1.0, 300., 1.0, 1.0]]
    g_p0s = [[1.0, 401., 0.1, 1.0, 403.2, 0.1],[1.0, 398.3, 0.1, 1.0, 401., 0.1],[1.0, 394.7, 0.1, 1.0, 397.3, 0.1]]

# gauss_exp_convo ranges
else:
    fit_type = 'gauss_exp_conv'
    n_ranges=[[331.5,346.5],[317.,331.],[290.5,307.5]]
    #g_ranges=[[400.,404.],[397.,402.5],[392.0,395.5]]
    g_ranges=[[395., 410.],[393., 408],[391, 402.]]
    n_p0s = [[1.0, 340., 1.0, 1.0],[1.0, 324., 1.0, 1.0],[1.0, 300., 1.0, 1.0]]
    g_p0s = [[1.0, 399., 0.1, 1.0, 401., 0.1, 1.0, 403.2, 0.1, 1.0, 406., 0.1],[1.0, 396., 0.1, 1.0, 398.3, 0.1, 1.0, 401., 0.1, 1.0, 403., 0.1],
             [1.0, 394., 0.1, 1.0, 394.7, 0.1, 1.0, 397.3, 0.1, 1.0, 400., 0.1]]

    g_bounds = (((-np.inf, 397, -np.inf, -np.inf, 400, -np.inf, -np.inf, 402, -np.inf, -np.inf, 405, -np.inf), 
                 (np.inf, 400, np.inf, np.inf, 402, np.inf, np.inf, 404, np.inf, np.inf, 407, np.inf)),
                ((-np.inf, 394, -np.inf, -np.inf, 397, -np.inf, -np.inf, 400, -np.inf, -np.inf, 402, -np.inf), 
                 (np.inf, 397, np.inf, np.inf, 399, np.inf, np.inf, 402, np.inf, np.inf, 404, np.inf)),
                ((-np.inf, 392, -np.inf, -np.inf, 394, -np.inf, -np.inf, 396, -np.inf, -np.inf, 398, -np.inf), 
                 (np.inf, 394.5, np.inf, np.inf, 396, np.inf, np.inf, 398, np.inf, np.inf, 402, np.inf)))


means_stds, max_n_time = [], []
for index,dist in enumerate(dists):
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = get_range(tof, 200, 450)
    n_tof = get_range(tof, n_ranges[index][0], n_ranges[index][1])
    g_tof = get_range(tof, g_ranges[index][0], g_ranges[index][1])
    
    # build hists  
    tof_hist, bin_centers = build_hist(tof, 1000)
    n_tof_hist, n_bin_centers = build_hist(n_tof, 1000)
    g_tof_hist, g_bin_centers = build_hist(g_tof, 1000)
    
    # neutrons
    p0 = n_p0s[index]
    if gauss == True:
        gmodel = lmfit.Model(gaussian)
        params = gmodel.make_params(a=1000,mu=p0[1],sigma=p0[2])
    
        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=True)
    
        coeff = (res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)
        n_hist_fit = gaussian(n_bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)
        n_full_hist_fit = gaussian(bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value)
        means_stds.append((coeff[1],coeff[2]))

    else:
        gmodel = lmfit.Model(gaus_exp_convo)
        params = gmodel.make_params(a=1000,mu=p0[1],sigma=p0[2],tau=p0[3])

        res = gmodel.fit(n_tof_hist, params, x=n_bin_centers, nan_policy='omit')
        print '\nFitting with LMFIT'
        print res.message
        print lmfit.fit_report(res,show_correl=True)

        coeff = (res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
        n_hist_fit = gaus_exp_convo(n_bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
        n_full_hist_fit = gaus_exp_convo(bin_centers, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value)
        means_stds.append((coeff[1]+coeff[3],coeff[2])) # gauss_exp_conv mean is mu + tau

    
    # get gauss exp conv maxes
    xvals = np.linspace(n_ranges[index][0], n_ranges[index][1], 1000)
    max_idx = np.argmax(gaus_exp_convo(xvals, res.params['a'].value, res.params['mu'].value, res.params['sigma'].value, res.params['tau'].value))
    max_n_time.append(xvals[max_idx])

    print "\nneutron gaussian fit:\n  mu = "+str(coeff[1])+" ns\n  sigma = "+str(coeff[2])

    # gammas
    p0 = g_p0s[index]
    bounds = g_bounds[index]
    #print p0
    #print bounds
    coeff, var_matrix = curve_fit(quad_gauss, g_bin_centers, g_tof_hist, p0=p0, bounds=bounds, max_nfev=10000)
    g_hist_fit = quad_gauss(g_bin_centers, *coeff)
    g_full_hist_fit = quad_gauss(bin_centers, *coeff)
    means_stds.append((coeff[1],coeff[2],coeff[4],coeff[5]))
    print_stats()
    
    # gamma plot
    fig1, ax1 = plt.subplots()
    plt.plot(g_bin_centers, g_tof_hist)
    plt.plot(g_bin_centers, g_hist_fit, '--') 
    plt.plot(g_bin_centers, gaussian(g_bin_centers, coeff[0], coeff[1], coeff[2]), '--')
    plt.plot(g_bin_centers, gaussian(g_bin_centers, coeff[3], coeff[4], coeff[5]), '--')
    plt.plot(g_bin_centers, gaussian(g_bin_centers, coeff[6], coeff[7], coeff[8]), '--')
    plt.plot(g_bin_centers, gaussian(g_bin_centers, coeff[9], coeff[10], coeff[11]), '--')

    plt.ylabel('counts')
    plt.xlabel('time (ns)')
    plt.text(0.7,0.75,'$\mu_1 = $'+str(round(means_stds[1+2*index][0],3))+'\n$\sigma_1 =$ '+str(round(means_stds[1+2*index][1],3))+'\n$\mu_2 =$ '+
             str(round(means_stds[1+2*index][2],3))+'\n$\sigma_2 =$ '+str(round(means_stds[1+2*index][3],3)), transform=ax1.transAxes)
    if plt_save:
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_g_fit.png',dpi=500)

    # neutron plot    
    fig2, ax2 = plt.subplots()
    plt.plot(n_bin_centers, n_tof_hist)
    plt.plot(n_bin_centers, n_hist_fit, label='convo')
    plt.ylabel('counts')
    plt.xlabel('time (ns)')
    plt.text(0.7,0.75,'$\mu =$ '+str(round(means_stds[0+2*index][0],3))+'\n$\sigma =$ '+str(round(means_stds[0+2*index][1],3)), transform=ax2.transAxes)
    if plt_save:    
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_n_fit.png',dpi=500)

    # full
    plt.figure()
    # scale
    n_full_hist_fit = [n*max(tof_hist)/max(n_full_hist_fit) for n in n_full_hist_fit]
    g_full_hist_fit = [g*max(tof_hist[len(tof_hist)*3/5:])/max(g_full_hist_fit) for g in g_full_hist_fit]
    plt.plot(bin_centers, tof_hist)
    #plt.plot(bin_centers, n_full_hist_fit, linewidth=2, linestyle='--')
    #plt.plot(bin_centers, g_full_hist_fit, linewidth=2, linestyle='--')
    plt.ylabel('Counts', fontsize=18)
    plt.xlabel('Time (ns)', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.xlim(325, 425)
    plt.tight_layout()
    plt.savefig(save_dir+fit_type+'_'+dist+'cm_tof_fits.pdf')
    if plt_save:
        plt.savefig(save_dir+fit_type+'_'+dist+'cm_tof_fits.png',dpi=500)

print '\n {:^8s} cm {:>8s} cm {:>8s} cm'.format(dists[0], dists[1], dists[2])
print '{:^8.2f} ns {:>8.2f} ns {:>8.2f} ns'.format(max_n_time[0], max_n_time[1], max_n_time[2])

plt.show()
if save_params:
    pickle.dump( means_stds, open( "peak_fit_params_11mev_"+fit_type+".p", "wb" ) )
    print '\nparams saved to peak_fit_params_11mev_'+fit_type+'.p'