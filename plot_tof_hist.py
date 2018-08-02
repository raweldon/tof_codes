import numpy as np
import matplotlib.pyplot as plt

#dists = ['180','256','363']
dists = ['179','276','369']

for dist in dists:
    tof_spec = np.load('dist_'+dist+'.npz')
    tof = tof_spec['data']
    tof = [x*4 for x in tof] # 1 clock cycle = 4 ns (250 MHz digitizer)
    tof = [i for i in tof if i<500 and i>200]

    plt.figure()
    plt.hist(tof,bins=1e4,histtype='step')
    plt.show()
