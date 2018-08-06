# tof_codes
tof codes for oct 2017 and feb 2018 runs

## scripts
tof_spec_fit.py - fits gaussians to neutron and gamma peaks   
tof_energy_calc.py - calculates energy and timing information

## additional info
8/2/18  
updated tof spectra on with daq2 codes   
had not been including all events (had cut by evno -- removed this)   
fixed multi peak issue for 11.3 MeV beam case   

8/3/18  
fixed issues with gauss_exp_conv (used as leading exp not trailing, found correct ranges)   
guassian and guass_exp_conv both look good   
update tof_energy_calc.py to account for gamma populations from havar and dcell stop 

8/6/18  
Working on 4 MeV branch
updated gamma fit to triple guassian with line for background -- looks good
REMEBER --> plot lmfit results with res.best_fit, intial with res.initial_fit 
   