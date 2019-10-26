# tof_codes
tof codes for oct 2017 and feb 2018 runs

## scripts
multi_exp_gauss_fit.py - used to calculate neutron peak time
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

8/9/18  
Working with multi neutron peaks    
tried using scipy optimize - not better than lmfit (still needed bounds)  
went back to lmfit - got gaussian model to work fairly well  

8/13/18  
Next steps:
-- gauss_exp_conv for gamma peaks  
-- better fit for 4 MeV neutron peaks  
-- error function for background of 11 MeV code  
-- update 11 MeV code to use lmfit  
-- find better functional fit for 11 MeV peak???  

8/15/2018  
Abandoned multi peak fits -- didn't really work  
Obtained neutron max peak location by eye, fits for gamma peaks
Results look great
   