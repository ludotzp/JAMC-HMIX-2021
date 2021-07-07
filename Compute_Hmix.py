import xarray as xr
import numpy as np
import os
from metpy import calc as mpcalc

#Directory and path to the sonde netcdf file
#The dropsonde dataset (George et al., 2021) can be downloaded here: https://doi.org/10.25326/221
#The radiosonde datasets (Stephan et al., 2020) can be downloaded here: https://doi.org/10.25326/137

input_dir = "~/Data/"
path_sondes = os.path.join(input_dir, "EUREC4A_JOANNE_Dropsonde-RD41_Level_3_v1.0.0.nc")

def load_sondes(path):
    
    #script to load the sondes and calculate relevant variables

    sondes = xr.open_dataset(path)
    
    sondes["q"].attrs['units'] = ''
    
    mixing_ratio = mpcalc.mixing_ratio_from_specific_humidity(sondes["q"])
    sondes["mixing_ratio"] = (["sonde_id", "alt"], mixing_ratio.magnitude)
        
    theta_v = mpcalc.virtual_potential_temperature(sondes["p"], sondes["ta"], sondes["mixing_ratio"])
    sondes["theta_v"] = (["sonde_id", "alt"], theta_v.magnitude)
    
    density = mpcalc.density(sondes["p"], sondes['ta'], sondes['mixing_ratio'])
    sondes["rho"] = (["sonde_id","alt"], density.magnitude)
           
    return sondes

def calculate_hmix(sondes):
    
    #Script to compute the height of the mixed layer in theta_v
       
    number_sondes = len(sondes.launch_time)
    hmix = np.zeros(number_sondes)
    thresh = 0.2 #in K, the threshold used to calculate the height of the mixed layer.
    number_min = 30 
    
    for i in range(number_sondes):

        sonde = sondes.isel(sonde_id=i).dropna(dim="alt",subset=["rho","theta_v"],\
                                                  how="any")     
        if (sonde.alt.min().values > 500): #if no measurements below 500 meters, we immediately drop the sonde
            sel_length = 0 
        else:      
            sel_length = len(sonde.alt.where(sonde.alt <= 500, drop=True))
            sonde = sonde.where(sonde.alt >= 100, drop=True)
       
        print(str(i)+"/"+str(number_sondes))

        #we keep only the sondes with at least 30 measurements of density and theta_v below 500 meters
        if (sel_length > number_min):
       
            var_thetav = 0
            numer_thetav = 0
            denom = 0
            thetav_mix = 0            
            k = 0 
           
            while(var_thetav < thresh):
                delta_z = sonde.alt[k+1]-sonde.alt[k]
                numer_thetav += 0.5*(sonde.rho[k+1]*sonde.theta_v[k+1] + sonde.rho[k]*sonde.theta_v[k])*delta_z
                denom += 0.5*(sonde.rho[k+1] + sonde.rho[k])*delta_z
                thetav_mix = numer_thetav/denom
                var_thetav = sonde.theta_v[k+1] - thetav_mix
                k += 1

            hmix[i] = sonde.alt.values[k]
            print(hmix[i])
            
            
        else:
            print(sel_length)
            print("sonde failed, hmix=0")
            hmix[i] = 0
    
    sondes["hmix"] = (("launch_time"), hmix)
    
    sondes = sondes.where(sondes.hmix > 0, drop=True)
       
    return sondes

sondes = load_sondes(path_sondes)
sondes = calculate_hmix(sondes)
#the modified netcdf file (with hmix) is saved in the input directory
sondes.to_netcdf(os.path.join(input_dir, "sondes_w_hmix.nc"))