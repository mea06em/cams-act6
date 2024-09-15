#!/usr/bin/python3

"""
author: felipe cifuentes (Adapted from John Douros)

Helper functions for comparing CAMS with TROPOMI
"""

import numpy as np
import xarray as xr


def CAMS_level_pressures(no2_cams, sp, a_inter, b_inter):
    
    #Dimensions will be use later
    nlev, nlat, nlon = no2_cams.shape
    
    pres_fl = np.zeros((nlev,nlat,nlon))
    pres_hl = np.zeros((nlev+1,nlat,nlon))

    # Calculate pressures at half and full levels
    for ilev in np.arange(nlev+1):

        # half-level pressures
        pres_hl[ilev,:,:] = a_inter[ilev] + b_inter[ilev] * sp

    for ilev in np.arange(nlev):
        # use half-level as & bs to produce full-level values
        a = (a_inter[ilev+1] + a_inter[ilev])/2.0
        b = (b_inter[ilev+1] + b_inter[ilev])/2.0
        pres_fl[ilev,:,:] = a + b * sp

    return pres_fl, pres_hl



def convert_o3_massmixing_to_molec_cm3(massmix, pressure, temperature):
    """ 
    Convert between kg/kg to molecule/m3

        * Pressures: Pa
        * Temperatures: K

        * Concentrations: μg/m3
        * Densities: molecules/cm3

    Ideal gas constant (sc.R) is in m3 Pa K−1 mol−1

    """ 
    import sys
    import scipy.constants as sc

    MWair = 28.97  # g/mol Molecular weight of dry air
    MWo3 = 48  # g/mol NO2
    ppm_to_molec_cm3 = 2.46e13  # Conversion factor for p = 1 atm and T = 298  (surface)

    volmix = massmix * MWair / MWo3

    volmix_ppm = volmix * 1.0e6

    dens = volmix_ppm * ( pressure / 101300.0 ) * ( 298.0 / temperature ) * ppm_to_molec_cm3
        
    return dens


def index_point(lat_point,lon_point,lat_list,lon_list):
    
    dif_lat = np.absolute(lat_list - lat_point)
    index_lat = dif_lat.argmin()
    
    dif_lon = np.absolute(lon_list - lon_point)
    index_lon = dif_lon.argmin()
    
    return index_lat, index_lon


def cams_to_tropomi_linear(no2_cams_profile, no2_tropomi_profile, pres_fl_cams, pres_fl_tropomi):

    nlev=len(no2_tropomi_profile)
    no2_cams_tropomi_grid=np.zeros((nlev))
    
    for x in range(nlev):
        
        dif_p = np.absolute(pres_fl_cams - pres_fl_tropomi[x])
        index_p1 = dif_p.argmin()
                            
        dif_p[index_p1] = dif_p[index_p1]*1e18
        index_p2 = dif_p.argmin()
                            
        m = (no2_cams_profile[index_p1] - no2_cams_profile[index_p2])/(pres_fl_cams[index_p1] - pres_fl_cams[index_p2])
        intercept = no2_cams_profile[index_p1] - m * pres_fl_cams[index_p1]
        
        no2_cams_tropomi_grid[x] = m * pres_fl_tropomi[x] + intercept
    
    return no2_cams_tropomi_grid