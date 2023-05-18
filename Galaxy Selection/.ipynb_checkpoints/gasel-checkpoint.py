# custom functions used frequently


import numpy as np
from astropy.coordinates import Distance
from astropy.cosmology import Planck18
from astropy import units as u
import astropy.constants as asc

def gal_ext(sky_coords):
    
    """
    Calculate the galactic extinction for an object, to be used in the flux correction
    Not yet implemented, and might never be implemented, unless needed
    """
    
    return 1


def app_mag(flux, mw_transmission=1):
    """
    Calculate the apparent magnitude from the flux as defined by https://www.legacysurvey.org/dr10/description/, including a mw transmission correction, if provided
    Input:
    - flux as numpy array in nanomaggies
    """
    
    return 22.5-2.5*np.log10(flux/mw_transmission)
    
    
    

def abs_mag(m, z, band, color_name, color):
    """
    Calculate the absolute magnitude, for now only uses m and z, no K-correction yet
    Input:
    - m: apparent magnitude
    - z: redshift 
    - band: band name, like "g"
    - color name: name of the color, e.g. "g - r"
    - color: value of the color
    """
    #return m - 5*np.log10(Distance(z=z, cosmology=Planck18)/u.Mpc*10**6)+5 - calc_kcor(band, z, color_name, color)
    return m - 5*np.log10(Distance(z=z, cosmology=Planck18)/u.Mpc*10**6)+5

# Calculating luminosities in different bands
# However, since we only look at the z-band we need to correct the solar luminosity
# Calculation from here: https://astronomy.stackexchange.com/questions/25126/how-to-calculate-luminosity-in-g-band-from-absolute-ab-magnitude-and-luminosity
lmbda_z = 920*10**(-9) #in  m
del_lambda_z = 160*10**(-9) #in m
del_v_z = (asc.c*lmbda_z/(del_lambda_z**2)).value # this is in /s now

lmbda_r = 640*10**(-9) #in  m 5600A-7200A -> 6400A median
del_lambda_r = 160*10**(-9) #in m
del_v_r = (asc.c*lmbda_r/(del_lambda_r**2)).value # this is in /s now


m_sun_z = -27.56 # in z-band: http://mips.as.arizona.edu/~cnaw/sun.html (DES filter)
m_sun_r = -27.12 # in r-band: http://mips.as.arizona.edu/~cnaw/sun.html (DES filter)

f_v_z = 10**((-48.6-m_sun_z)/2.5) # in erg/(cm^2 s Hz)
f_v_r = 10**((-48.6-m_sun_r)/2.5) # in erg/(cm^2 s Hz)

dist_sun = 1.496*10**13 # in cm

L_sun_z = f_v_z*del_v_z*4*np.pi*dist_sun**2 *10**(-7) # this is  in Watts now
L_sun_r = f_v_r*del_v_r*4*np.pi*dist_sun**2 *10**(-7)
# the value is roughly half of the full bolometric value... not sure if this makes sense

def lum_z(M):
    """
    Get the luminosity in the z-band
    Input:
    - absolute magnitude M in this band
    """
    return L_sun_z*10**(-0.4*M)*u.W

def lum_r(M):
    """
    Get the luminosity in the r-band
    Input:
    - absolute magnitude M in this band
    """
    return L_sun_r*10**(-0.4*M)*u.W

def lum_wiggle(D_gal, m):
    """
    Calculate a luminositiy proxy from the magnitude
    PArams:
    - D_gal: Distance to the galaxy, measured from redshift
    - m: magnitude (for us in z-band)
    
    returns: D_gal**2*10**(-0.4*m)
    """
    return D_gal**2*10**(-0.4*m)


def V_max(omega_s, z_min, z_max_lum):
    """
    Description: Calculates the maximum Volume in which the source could have been detected in: corrects for the so-called "Malmquist-bias"
    (faint objects, which usually also means low mass objects, will only be covered in a survey within a smaller volume than bright and high-mass
    objects)
    
    Params:
    omega_s: surface area covered by the complete data
    z_min: lower redshift limit of sample
    z_max_lum: maximum redshift determined for for object based on luminosity and its completeness 
    (i.e. maximum redshift at which the source with luminosity l would no longer be part of the sample)
    
    """
    z_max = np.minimum(z_max_data, z_max_lum)
    return 4/3*np.pi * omega_s/omega_sky * (Distance(z=z_max, cosmology=Planck18)**3-Distance(z=z_min, cosmology=Planck18)**3)


def lum_lim(lum, m_lim, m):
    """
    Calculate the luminosity a source would have if its magnitude was equal to the DESI magnitude limit (have to research the value)
    
    Params:
    - lum: luminosity of an object
    - m_lim: magnitude limit of the survey
    - m: apparent magnitude of an object
    
    """
    return lum*10**(-0.4*(m_lim-m))



def calc_kcor(filter_name, redshift, colour_name, colour_value):
    """
    K-corrections calculator in Python. See http://kcor.sai.msu.ru for the 
    reference. Available filter-colour combinations must be present in the 
    `coeff` dictionary keys.

    @type   filter_name: string    
    @param  filter_name: Name of the filter to calculate K-correction for, e.g. 
                         'u', 'g', 'r' for some of the SDSS filters, or 'J2', 
                         'H2', 'Ks2' for 2MASS filters (must be present in 
                         `coeff` dictionary)
    @type      redshift: float    
    @param     redshift: Redshift of a galaxy, should be between 0.0 and 0.5 (no
                         check is made, however)
    @type   colour_name: string    
    @param  colour_name: Human name of the colour, e.g. 'u - g', 'g - r', 
                         'V - Rc', 'J2 - Ks2' (must be present in `coeff` dictionary)
    @type  colour_value: float    
    @param colour_value: Value of the galaxy's colour, specified in colour_name    
    @rtype:              float
    @return:             K-correction in specified filter for given redshift and 
                         colour
    @version:            2012
    @author:             Chilingarian, I., Melchior. A.-L., and Zolotukhin, I.
    @license:            Simplified BSD license, see http://kcor.sai.msu.ru/license.txt

    Usage example:
    
        >>> calc_kcor('g', 0.2, 'g - r', 1.1)
        0.5209713975999992
        >>> calc_kcor('Ic', 0.4, 'V - Ic', 2.0)
        0.310069919999993
        >>> calc_kcor('H', 0.5, 'H - K', 0.1)
        -0.14983142499999502
        
    """
    coeff = {

        'r_gr': [
            [0,0,0,0],
            [1.83285,-2.71446,4.97336,-3.66864],
            [-19.7595,10.5033,18.8196,6.07785],
            [33.6059,-120.713,-49.299,0],
            [144.371,216.453,0,0],
            [-295.39,0,0,0],
        ], 
        'g_gz': [
            [0,0,0,0],
            [2.37454,-4.39943,7.29383,-2.90691],
            [-28.7217,-20.7783,18.3055,5.04468],
            [220.097,-81.883,-55.8349,0],
            [-290.86,253.677,0,0],
            [-73.5316,0,0,0],
        ],
        'z_gz': [
            [0,0,0,0],
            [0.30146,-0.623614,1.40008,-0.534053],
            [-10.9584,-4.515,2.17456,0.913877],
            [66.0541,4.18323,-8.42098,0],
            [-169.494,14.5628,0,0],
            [144.021,0,0,0],
        ]

    }

    c = coeff[filter_name + '_' + colour_name.replace(' - ', '')]
    kcor = 0.0

    for x, a in enumerate(c):
        for y, b in enumerate(c[x]):
            kcor += c[x][y] * redshift**x * colour_value**y

    return kcor