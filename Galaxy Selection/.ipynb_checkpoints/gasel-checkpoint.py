# custom functions used frequently


import numpy as np
from astropy.coordinates import Distance
from astropy.cosmology import Planck18
from astropy import units as u
import astropy.constants as asc
from astropy.table import Table, hstack, vstack
from astropy import table


import psycopg2


omega_sky = 41253 #

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
L_sun_bol = 3.828*10**26

L_sun_bands = {"r": L_sun_r*u.W, "z": L_sun_z*u.W, "bol": L_sun_bol*u.W}

std_names = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'LASTNIGHT', 'Z', 'ZERR', 'ZWARN', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'SPECTYPE', 'BGS_TARGET', 'EBV', 'SERSIC']

ext_coeffs =  [3.995, 3.214, 2.165, 1.592, 1.211, 1.064] # for the DECam 𝑢, 𝑔, 𝑟, 𝑖, 𝑧, 𝑌 

def connect_to_time_db():
    """
    connect to the desi daily SQL database
    
    Returns
    -------
        cursor: cursor object
            cursor for the database
    """
    try:
        db = psycopg2.connect(host='decatdb.lbl.gov', database='desidb', user='desi', password = "5kFibers!", port="5432")
        cursor = db.cursor()
        return cursor
    except (Exception, psycopg2.Error) as error:
        print("FAILED TO ESTABLISH CONNECTION TO DATABASE")
        print("------------------------------------------")
        print("Error log: ")
        print(error)
        
def make_table(rows, names = std_names):
    """
    Creates an astropy data table with only good redshifts and only unique entries (takes the first)
    
    Parameters
    ----------
        rows: row object from cursor.fetchall() for example
            contains the actual data from the sql database 
        names: list
            list of names for the columns (e.g. ["TARGET_RA", "TARGET_DEC"])
            
    Returns
    -------
        table: astropy table
            returns an astropy table with the specified column names and only good redshifts; also ensures only unique entries exist and can return a table with no data, only column names
    """
    if rows:
        data = Table(list(map(list, zip(*rows))),
                             names=names)
        
        if len(data) > 0:
            data = data[data['ZWARN']==0]
            data = data[data['Z']>=0]
            data = table.unique(data, keys = "TARGETID")
        return data
    else:
        return Table(names = names)

def make_list(rows, names = std_names):
    """
    Creates a numpy list from the sql rows
    
    Parameters
    ----------
        rows: row object from cursor.fetchall() for example
            contains the actual data from the sql database
            
    Returns
    -------
        data: python list
            returns a python list without any data clean up
    """
    
    if rows:
        data = list(map(list, zip(*rows)))
        return data
    

def get_data_rad_search(ra, dec, radius, cursor, query=None):
    
    """
    Queries data from the database; a standard query is implemented, but can be altered
    Parameters
    ----------
        ra: int/float
            right ascension of the center for the circular search
        dec: int/float
            declination of the center of the circular search
        radius: int/float
            search radius
        cursors: db.cursor()
            cursors object from the database
        query (optional): str
            the query string (i.e. what data to retrieve from the database
    
    Returns
    -------
        rows: row object from database
            contains the actual information. use make_table() to convert to astropy table
            
    """
    redux = "daily"
    if query == None:
        query = 'SELECT f.targetid, f.target_ra, f.target_dec, c.night, r.z, r.zerr, r.zwarn, f.flux_r, f.flux_g, f.flux_z, r.spectype, f.bgs_target, f.ebv, f.sersic\n' \
                    f'FROM {redux}.tiles_fibermap f\n' \
                    f'INNER JOIN {redux}.cumulative_tiles c ON f.cumultile_id=c.id\n' \
                    f'INNER JOIN {redux}.tiles_redshifts r ON r.cumultile_id=c.id AND r.targetid=f.targetid\n' \
                    f'WHERE q3c_radial_query( f.target_ra, f.target_dec, {ra}, {dec}, {radius});'
    cursor.execute(query)
    rows = cursor.fetchall()
    return rows
       
def db_doall(ra, dec, radius, names = std_names, query=None):
    
    """
    Do all steps to retrieve data from the database. A standard query is provided, but can be changed
    
    Parameters
    ----------
        ra: int/float
            right ascension of the center for the circular search
        dec: int/float
            declination of the center of the circular search
        radius: int/float
            search radius
        names (optional): list
            list of names for the columns (e.g. ["TARGET_RA", "TARGET_DEC"])
        query (optional): str
            the query string (i.e. what data to retrieve from the database
    Returns
    -------
        table: astropy table
            returns the data table with only good redshift objects and unique entries
    """
    
    cursor = connect_to_time_db()
    rows = get_data_rad_search(ra, dec, radius, cursor, query)
    table = make_table(rows, names)
    return table
    
    
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

def M(L, band = "bol"):
    
    """
    Calculates the absolute magnitude from the luminosity
    One can specify a band or just use the bolometric value without specifiyng a band
    Input:
    - L: luminosity in Watts (u.W)
    - band: "g", "r", "z" for example or specify nothing for bolometric
    """
    
    L_sun = L_sun_bands[band]
    
    return -2.5*np.log10(L/(L_sun))


def lum(M, band = "bol"):
    """
    Get the luminosity in watts
    Input:
    - absolute magnitude M in this band
    """
    L_sun = L_sun_bands[band]
    return L_sun*10**(-0.4*M)

def lum_wiggle(D_gal, m):
    """
    Calculate a luminositiy proxy from the magnitude
    PArams:
    - D_gal: Distance to the galaxy, measured from redshift
    - m: magnitude (for us in z-band)
    
    returns: D_gal**2*10**(-0.4*m)
    """
    return D_gal**2*10**(-0.4*m)


def V_max(omega_s, z_min, z_max_data, z_max_lum):
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