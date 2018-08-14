''' Caculates TOF and beam energy using params file generated with tof_spec_fit.py
    Assumes gamma flash is from havar foil (distance is from deuterium cell through point of interaction in EJ-309 2x2 det)
'''


#!/usr/bin/env python
import numpy as np
import pickle

def full_calc(dist,tof):
    # calculated using E = 0.5mv^2 (constants = 0.5227)
    return mn/2./c**2*(dist/tof)**2

def calc_dist(tof):
    # calculate distance assuming energy is known
    return tof*np.sqrt(E_n/(mn/2./c**2))

def calc_tof(dist):
    return dist/np.sqrt(E_n/(mn/2./c**2))

def rel_calc(dist,tof,m0):
    v = dist/tof
    gamma = 1/np.sqrt(1-(v/c)**2)
    ke = (gamma-1)*m0
    return ke

def rel_dist(tof):
    return tof*c*np.sqrt(1.-(1./(1.+(E_n/mn)))**2)

def rel_tof(dist):
    return dist/(c*np.sqrt(1.-(1./(1.+(E_n/mn)))**2))

def nonrel_resolution(dist, tof, del_tof):
    return 2*E_n*np.sqrt((del_dist/dist)**2 + (del_tof/tof)**2)

def rel_resolution(dist,tof,del_tof):
    beta = dist/tof/c
    res = (E_n + mn)/E_n * beta**2/(1-beta**2) * np.sqrt((del_dist/dist)**2+(del_tof/tof)**2)
    return res     

def fom(dist,tof):
    # want a value ~23 ns/m (tsoulifinidis)
    return 1/(np.sqrt(1 - (mn/(E_n+mn))**2)*c)
    
def nonrel_vel(T,m):
    return np.sqrt(2.*T/m*c**2)

def rel_vel(T,m):
    return c*np.sqrt(1. - (1./(1.+T/m))**2) #cm/ns

peak_params = pickle.load( open('peak_fit_params_11mev_gauss.p','r'))

tof_n = [peak_params[0][0],peak_params[2][0],peak_params[4][0]]
print tof_n

# from gcrich -- gamma peaks (1) unknown, (2) collimator, (3) havar, (4) d cell backing
tof_g1 = [peak_params[1][0],peak_params[3][0],peak_params[5][2]] # havar
tof_g2 = [peak_params[1][2],peak_params[3][2],peak_params[5][0]] # d cell backing
sigma_g1 = [peak_params[1][1],peak_params[3][1],peak_params[5][3]] # sigma of gamma peak
sigma_g2 = [peak_params[1][3],peak_params[3][3],peak_params[5][1]] # sigma of gamma peak
print tof_g1

# constants
dist = [180.6,255.6,362.8] # dist from det to end of beam line
dcel_length = 2.8575 # cm
det_size = 5.08 # cm, avg interaction depth in 0 deg det
dcel_to_beamend = 145.2 # 146.58 measured by Ron Malone, 145.2 from our measurement (2-18)
c = 29.9792 # cm/ns
mn = 939.552 # MeV neutron mass  
E_n = 11.325 # MeV
md = 1875.6 # MeV
E_d_start = 8.26 # E after havar (8.55 MeV before foil, 8.26 after, 8.209 at cell center, 0.103 MeV loss in cell)
E_d_14 = 8.26 - 0.103*0.25
E_d_34 = 8.26 - 0.103*0.75  # 3/4 through d cell

# uncerts
del_dist = np.sqrt((dcel_length/2.)**2 + (5.08/2.)**2) # 5.08cm len of det
bunch_width = [2.355*s for s in sigma_g1] # fwhm of gamma peak gives estimation of bunch width
det_traverse = 5.08/rel_vel(E_n,mn)
del_tof_rel = [np.sqrt(det_traverse**2 + b**2) for b in bunch_width]
del_tof_nonrel = [np.sqrt((5.08/nonrel_vel(E_n,mn))**2 + b**2) for b in bunch_width]

havar_peak = True

if havar_peak:
    dist_to_dcell = [dcel_to_beamend + dcel_length/2. + det_size/2. + d for d in dist] # for havar gamma start
    #     gamma peak            gamma travel time from havar to detection      n peak  deuteron travel time to center of havar
    tof_rel = [tof_g1[i] + (dist_to_dcell[i] + dcel_length/2.)/c - n - dcel_length/2./rel_vel(E_d_14,md) for i,n in enumerate(tof_n)] # havar gamma start
    tof_nonrel = [tof_g1[i] + (dist_to_dcell[i] + dcel_length/2.)/c - n - dcel_length/2./nonrel_vel(E_d_14,md) for i,n in enumerate(tof_n)]

else:
    dist_to_dcell = [dcel_to_beamend + det_size/2. + d for d in dist] # for dcell stop gamma start
    #     gamma peak   gamma travel time from dcell stop to det    n peak  n travel time from dcell center to stop
    tof_rel = [tof_g2[i]+(dist_to_dcell[i])/c + dcel_length/2.*rel_vel(E_d_34,md) - n     for i,n in enumerate(tof_n)] # dcell stop gamma start
    tof_nonrel = [tof_g2[i]+(dist_to_dcell[i])/c + dcel_length/2.*nonrel_vel(E_d_34,md) - n     for i,n in enumerate(tof_n)]

dists = ['180 cm', '256 cm', '363 cm']
for i,t in enumerate(tof_rel):
    print '--------------------------------------'
    print dists[i],'\n'
    print 'non rel En      =', full_calc(dist_to_dcell[i],tof_nonrel[i]), ' MeV'
    print 'measured dist   =', dist_to_dcell[i]
    print 'calculated dist =', calc_dist(tof_nonrel[i]), ' cm'
    print 'difference      =', calc_dist(tof_nonrel[i]) - dist_to_dcell[i], ' cm'
    print 'calculated tof  =', calc_tof(dist_to_dcell[i]), ' ns'
    print 'measured tof    =', tof_nonrel[i]
    print 'difference      =', calc_tof(dist_to_dcell[i]) - tof_nonrel[i], ' ns'
    print 'resolution      =', nonrel_resolution(dist_to_dcell[i],tof_nonrel[i],del_tof_nonrel[i]), ' MeV\n'
    print 'rel En          =', rel_calc(dist_to_dcell[i],t,mn), ' MeV'
    print 'calculated dist =', rel_dist(t), ' cm'
    print 'measured dist   =', dist_to_dcell[i], ' cm'
    print 'difference      =', rel_dist(t) - dist_to_dcell[i], ' cm'
    print 'calculated tof  =', rel_tof(dist_to_dcell[i]), ' ns'
    print 'measured tof    =', t, ' ns'
    print 'tof difference  =', rel_tof(dist_to_dcell[i]) - t, ' ns'
    print 'resolution      =', rel_resolution(dist_to_dcell[i],t,del_tof_rel[i]
           ), '\n                =', rel_resolution(dist_to_dcell[i],t,del_tof_rel[i])*E_n, ' MeV'
    print 'FOM             =', fom(dist_to_dcell[i],t),'\n'

