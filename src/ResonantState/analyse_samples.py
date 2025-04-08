import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import constants as c
import warnings
warnings.filterwarnings('ignore')

def convert_units(planet_values, param, mstar, rstar, units):
    """
    Converts planetary values from stellar units to the specified system.

    Parameters:
    -----------
    planet_values : array-like
        Values to convert (in stellar units).
    param : str
        Parameter type ('mass', 'radius', or 'density').
    mstar : array-like
        Stellar mass in solar masses.
    rstar : array-like
        Stellar radius in solar radii.
    units : str
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    array-like
        Converted parameter values in the desired units.
    """

    m_star_kg, r_star_m = mstar * c.M_sun, rstar * c.R_sun
    rho_star_SI = (m_star_kg / r_star_m**3)

    if param in ['mass','radius','density'] and units!='star':
        unit_map = {
            'mass' : {
                'star_to_SI': m_star_kg,
                'earth': c.M_earth,
                'jup' : c.M_jup,
                'sun': c.M_sun,
                'SI': 1.
            },
            'radius': {
                'star_to_SI': r_star_m,
                'earth': c.R_earth,
                'jup' : c.R_jup,
                'sun': c.R_sun,
                'SI': 1.            
            },
            'density':{
                'star_to_SI': rho_star_SI,
                'earth': c.M_earth/(c.R_earth**3),
                'jup' : c.M_jup/(c.R_jup**3),
                'sun': c.M_sun/(c.R_sun**3),
                'SI': 1000.       
            }
        }
        planet_values_SI = planet_values * unit_map.get(param, {}).get('star_to_SI', 1)
        return planet_values_SI / unit_map[param].get(units, 1)
    else:
        return planet_values


def get_parameter_key(param, p):
    par_keys = {
                'mass': f'mass_planet_star_ratio_{p}',
                'ecc': [f'k_{p}', f'h_{p}'],
                'radius': f'radius_planet_star_ratio_{p}',
                'period': f'period_days_{p}',
                'varpi': [f'h_{p}', f'k_{p}'],
                'k': f'k_{p}',
                'h': f'h_{p}',
                'omega': f'longitude_of_ascending_node_deg_{p}',
                'lambda': f'mean_longitude_deg_{p}',
                'incl': f'inclination_deg_{p}',
                'density': [f'mass_planet_star_ratio_{p}', f'radius_planet_star_ratio_{p}']
                }
    return par_keys[param]

def get_labels(param, units):
    par_labels = {
        'mass': fr'$M_{{\mathrm{{planet}}}} \; (M_{{\mathrm{{{units}}}}})$',
        'ecc': r'$e$',
        'radius': fr'$R_{{\mathrm{{planet}}}} \; (R_{{\mathrm{{{units}}}}})$',
        'period': r'$P_{\mathrm{planet}} \; (d)$',
        'varpi': r'$\varpi \; (deg)$',
        'density': r'$g/cm^3$' if units=='SI' else fr'$\rho_{{\mathrm{{planet}}}} \; (\rho_{{\mathrm{{{units}}}}})$',
        'omega': r'$\Omega \; (deg)$',
        'lambda': r'$\lambda \; (deg)$',
        'incl': 'i (deg)',
        'k': r'$k = ecos(\varpi)$',
        'h': r'$h = esin(\varpi)$'
    }
    return par_labels[param]


def get_samples(df, param,  p, units):
    """
    Extracts and converts samples for a specific parameter and planet.

    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe containing posterior samples.
    param : str
        Name of parameter to be extracted.
    p : int
        Planet index.
    units : str
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    numpy.ndarray
        Converted samples for the selected parameter and planet.
    """

    param_key = get_parameter_key(param,p)

    if param in ['ecc', 'varpi', 'density']:
        par1, par2 = df[param_key[0]].values, df[param_key[1]].values
        if param=='ecc':
            samples = np.sqrt(par1**2 + par2**2)
        if param=='varpi':
            samples = np.arctan2(par1,par2) * (180 / np.pi)  
        if param=='density':
            # Avoid division by zero or near-zero values
            valid = par2 > 1e-6  
            if not np.any(valid):   # if all values are invalid
                return np.array([])  
            samples = par1[valid] / (par2[valid]**3)
            mstar = df['mass_star_m_sun'].values[valid]
            rstar = df['radius_star_r_sun'].values[valid]
            return convert_units(samples, param, mstar, rstar, units)
    else:
        samples = df[param_key].values

    mstar = df['mass_star_m_sun'].values
    rstar = df['radius_star_r_sun'].values

    samples_in_units = convert_units(samples,param, mstar, rstar, units)
    return samples_in_units


def get_nb_planets(data):
    nb_planets = 0
    for df_dict in data:
        n_planets = len(df_dict["planets_list"])
        if n_planets > nb_planets:
            nb_planets = n_planets
    return nb_planets
        

def plot_histograms(dict_list, param, units='star'):
    """
    Plots histograms of a given parameter for each planet in each unique analysis.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    param : str
        Parameter to plot.
    units : str, optional (default is 'star')
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI').
        Applies only to parameters 'mass', 'radius', or 'density'.
    """

    nb_planets = get_nb_planets(dict_list)
    fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
    axs = np.atleast_1d(axes)

    # Loop through each analysis and plot the histograms
    for df_dict in dict_list:
        key = df_dict['sample_name']
        df = df_dict['sample']

        for p in range(nb_planets):
            if p < len(df_dict['planets_list']):
                x = get_samples(df, param,  p, units)
                axs[p].hist(x, bins=50, alpha=0.5, label=key)
                axs[p].set_xlabel(get_labels(param, units))
                axs[p].legend()
                axs[p].set_title(f'planet {p}')
    plt.tight_layout()
    plt.show()


def plot_samples(dict_list, x_param, y_param, units='star'):
    """
    Plots scatter plots of two parameters for each planet in each unique analysis.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    x_param : str
        Parameter for the x-axis.
    y_param : str
        Parameter for the y-axis.
    units : str, optional (default is 'star')
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI'). 
        Applies only to parameters 'mass', 'radius', or 'density'.
    """

    nb_planets = get_nb_planets(dict_list)
    fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
    axs = np.atleast_1d(axes)

    # Loop through each analysis and plot the samples
    for df_dict in dict_list:
        key = df_dict['sample_name']
        df = df_dict['sample']

        for p in range(nb_planets):
            if p < len(df_dict['planets_list']):
                axs[p].set_xlabel(get_labels(x_param, units))
                axs[p].set_ylabel(get_labels(y_param, units))
                axs[p].legend()
                axs[p].set_title(f'planet {p}')

                x = get_samples(df, x_param,  p, units)
                y = get_samples(df, y_param,  p, units)
                if (x.size == 0) or (y.size == 0):
                    continue
                axs[p].scatter(x, y, alpha=0.1, label=key)

    plt.tight_layout()
    plt.show()


def compare_period_ratios(dict_list):
    """
    Plots histograms of period ratios (P_{i+1} / P_i) between adjacent planets.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    """

    nb_planets = get_nb_planets(dict_list)
    fig, axes = plt.subplots(nb_planets-1, 1, figsize=(6, 5 * (nb_planets-1)))
    axs = np.atleast_1d(axes)

    # Loop through each analysis and plot the period ratios
    for df_dict in dict_list:
        key = df_dict['sample_name']
        df = df_dict['sample']

        for p in range(nb_planets - 1):
            if p + 1 < len(df_dict['planets_list']):
                val_p = get_samples(df, 'period',  p, units='star')
                val_p1 = get_samples(df, 'period',  p+1, units='star')
                axs[p].hist(val_p1 / val_p, bins=50, alpha=0.5, label=key)
                axs[p].set_xlabel(f'P_{p+1} / P_{p}')
                axs[p].legend()
    plt.tight_layout()
    plt.show()


def plot_consecutive_planets(dict_list, param, units='star'):
    """
    Plots scatter plots comparing a parameter between consecutive planets.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    param : str
        Parameter to compare.
    units : str, optional (default is 'star')
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI'). 
        Applies only to parameters 'mass', 'radius', or 'density'.
    """
    nb_planets = get_nb_planets(dict_list)
    fig, axes = plt.subplots(nb_planets-1, 1, figsize=(6, 5 * (nb_planets-1)))
    axs = np.atleast_1d(axes)

    # Loop through each analysis and plot the samples of consecutive planets
    for df_dict in dict_list:
        key = df_dict['sample_name']
        df = df_dict['sample']

        for p in range(nb_planets - 1):
            if p + 1 < len(df_dict['planets_list']):
                axs[p].set_xlabel(f'{get_labels(param, units)} (planet {p})')
                axs[p].set_ylabel(f'{get_labels(param, units)} (planet {p+1})')
                axs[p].set_title(f'{param}_{p+1} vs {param}_{p}')
                axs[p].legend()

                x = get_samples(df, param,  p, units)
                y = get_samples(df, param,  p+1, units)
                if (x.size == 0) or (y.size == 0):
                    continue
                axs[p].scatter(x, y, alpha=0.1, label=key)

    plt.tight_layout()
    plt.show()
