from scipy.optimize import curve_fit
import math
import copy
import os
import netCDF4
import pandas as pd
from datetime import timedelta
from fnmatch import fnmatch
import numpy as np
import itertools
from scipy import signal
import qc_utils

# constants from Nimbus
mol_wgt_dry_air = 28.9637 # kg/kmol
R0 = 8314.462618 # J/kmol/K
Rd = R0/mol_wgt_dry_air
Cpd = 7.0/2.0*Rd
Cvd = 5.0/2.0*Rd

def calc_mach_dry(q: np.array, ps: np.array):
    # a function to calculate mach given a dynamic (q) and static (ps) pressure
    # Sometimes q is small and negative. When it's negative, mach is imaginary. So, only calc
    # mach when q is >= 0.
    mach_dry = np.zeros(len(q))
    q_pos_inds = q >= 0
    mach_dry[q_pos_inds] = (5*(((ps[q_pos_inds] + q[q_pos_inds])/ps[q_pos_inds])**(Rd/Cpd) - 1))**0.5
    return mach_dry

def calc_akrd(x: np.array, a: float, b: float, c: float, q: np.array = None):
     # x: a 2 by N array, where x[0,:] = ADIFR/QCF and x[1,:] = dry mach number
     # a, b, c: fit coefficients
     ratio = x[0,:]
     mach_dry = x[1,:]
     akrd = np.zeros(len(ratio))
     if q is not None:
         valid = q > 5.5 # true where dynamic pressure is greater than 5.5 hPa, like in nimbus
     else:
         valid = np.array(np.ones(len(ratio)), dtype='bool')
     akrd[valid] = a + ratio[valid]*(b + c*mach_dry[valid])
     return akrd
    

def calc_winds(data_df: pd.DataFrame, aoa_coefs: np.array, aos_coefs: np.array, name_append: str = '_test') -> pd.DataFrame:

    # first, calculate AoA and AoS with given coefficients
    a_ratio = data_df['ADIFR']/data_df['QCF']
    b_ratio = data_df['BDIFR']/data_df['QCF']
    mach = qc_utils.calc_mach(data_df['QCF'], data_df['PSF'])
    akrd_new = qc_utils.fit_func(np.array([a_ratio,mach]), *aoa_coefs)
    sslip_new = qc_utils.simple_fit_func(b_ratio, *aos_coefs)
    # add them to the data frame
    data_df['AKRD'+name_append] = akrd_new
    data_df['SSLIP'+name_append] = sslip_new


