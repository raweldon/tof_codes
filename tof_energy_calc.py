#!/usr/bin/env python
import numpy as np
import pickle

def full_calc(dist,tof):
    # calculated using E = 0.5mv^2 (constants = 0.5227)
    return n_mass/2*u_to_mev/c**2*(dist/tof)**2

def calc_dist(tof):
    # calculate distance assuming energy is known
    return tof*np.sqrt(E_n/(n_mass/2*u_to_mev/c**2))

def calc_tof(dist):
    return dist/np.sqrt(E_n/(n_mass/2*u_to_mev/c**2))

def rel_calc(dist,tof,m0):
    v = dist/tof
    #print v
    gamma = 1/np.sqrt(1-(v/c)**2)
    ke = (gamma-1)*m0
    return ke

def rel_dist(tof):
    return tof*c*np.sqrt(1.-(1./(1.+(E_n/mn)))**2)

def rel_tof(dist):
    return dist/(c*np.sqrt(1.-(1./(1.+(E_n/mn)))**2))
    
def resolution(dist,tof,del_tof):
    beta = dist/tof/c
    res = (E_n + mn)/E_n * beta**2/(1-beta**2) * np.sqrt((del_dist/dist)**2+(del_tof/tof)**2)
    return res

def fom(dist,tof):
    # want a value ~23 ns/m (tsoulifinidis)
    return 1/(np.sqrt(1 - (mn/(E_n+mn))**2)*c)
    
def rel_vel(T,m):
    return np.sqrt(c*(1.-(1./(1.+T/m))**2))

peak_params = pickle.load( open('peak_fit_params.p','rb'))

tof_n = [peak_params[0][0],peak_params[2][0],peak_params[4][0]]
tof_g1 = [peak_params[1][2],peak_params[3][0],peak_params[5][2]] # first gamma peak (in time), probably from colimator
tof_g2 = [peak_params[1][0],peak_params[3][2],peak_params[5][0]] # second gamma peak (assuming this is from d cell interactions)

dist = [175.8,234.9,284.3] # dist from det to end of beam line
dcel_length = 2.8575 # cm
ej309_interaction = 2. # cm, avg interaction depth in 0 deg det
dcel_to_beamend = 146.58 # measured by Ron Malone, 145.2 from our measurement
dist_to_dcell = [dcel_to_beamend+dcel_length/2+ej309_interaction+d for d in dist] 
n_mass = 1.008664 # u
u_to_mev = 931.4941 # u = 931 MeV/c^2
c = 29.9792 # cm/ns
mn = 939.552 # MeV neutron mass  
E_n = 11.325 # MeV
md = 1875.6 # MeV
E_d = 8.3 # E after havar
d_vel = rel_vel(E_d,md)

#del_dist = np.sqrt(dcel_length**2+5.08**2) # 5.08cm len of det
#del_tof = np.sqrt(3**2+2**2) # 3ns for spread in n peak, 2ns for spread in gamma peak
del_dist = dcel_length/2.+5.08/2. # 5.08cm len of det
#del_tof = [peak_params[0][1]+peak_params[1][1],
#           peak_params[2][1]+peak_params[3][3],
#           peak_params[4][1]+peak_params[5][1]]
del_tof = [0,0,0] # uncert should be sigma/sqrt(N) -> uncert in the mean is very small (<0.01); most of uncert is due to the uncert in distance
#tof = [tof_g2[i]+dist_to_dcell[i]/c - n for i,n in enumerate(tof_n)] #gamma peak time + gamma travel time - neutron peak time 
#print tof
#     gamma peak            gamma travel time from havar to detection       n peak  deuteron travel time to center of havar
tof = [tof_g2[i]+(dcel_to_beamend+dcel_length+ej309_interaction+dist[i])/c - n - d_vel*dcel_length/2 for i,n in enumerate(tof_n)]
#print tof

dists = ['175 cm', '235 cm', '284 cm']
for i,t in enumerate(tof):
    print '--------------------------------------'
    print dists[i]
    print 'non rel En      =', full_calc(dist_to_dcell[i],t), ' MeV'
    print 'calculated dist =', calc_dist(t), ' cm'
    print 'difference      =', calc_dist(t) - dist_to_dcell[i], ' cm'
    print 'calculated tof  =', calc_tof(dist_to_dcell[i]), ' ns'
    print 'difference      =', calc_tof(dist_to_dcell[i]) - t, ' ns\n'
    print 'rel En          =', rel_calc(dist_to_dcell[i],t,mn), ' MeV'
    print 'calculated dist =', rel_dist(t), ' cm'
    print 'measured dist   =', dist_to_dcell[i], ' cm'
    print 'difference      =', rel_dist(t) - dist_to_dcell[i], ' cm'
    print 'calculated tof  =', rel_tof(dist_to_dcell[i]), ' ns'
    print 'measured tof    =', t, ' ns'
    print 'tof difference  =', rel_tof(dist_to_dcell[i]) - t, ' ns'
    print 'resolution      =', resolution(dist_to_dcell[i],t,del_tof[i]
           ), '\n                =', resolution(dist_to_dcell[i],t,del_tof[i])*E_n, ' MeV'
    print 'FOM             =', fom(dist_to_dcell[i],t),'\n'

