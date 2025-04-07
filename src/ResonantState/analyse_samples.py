import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import constants as c

def convert_units(planet_values, param, mstar, rstar, units):
    """
    Converts planet parameter values (mass, radius, or density) from stellar units 
    to the specified target unit system.

    Parameters:
    -----------
    planet_values : array-like
        Array of values in stellar units (e.g., mass or radius ratios).
    param : str
        The parameter type to convert ('mass', 'radius', or 'density').
    mstar : float or array-like
        Mass of the host star in solar masses.
    rstar : float or array-like
        Radius of the host star in solar radii.
    units : str
        Target unit system ('earth', 'jup', 'sun', 'SI', or 'star').

    Returns:
    --------
    converted_values : array-like
        Values converted into the specified unit system.
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
        'incl': 'i (deg)'
    }
    return par_labels[param]


def extract_samples(df, param,  p, units):
    """
    Extracts posterior samples from input dataframe for a specified parameter and planet index.
    Converts the extracted values into the specified units.

    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe containing the posterior samples from one analysis.
    param : str
        Parameter to extract.
    p : int
        Planet index starting from 0.
    units : str
        Unit system to convert the mass/radius/density samples ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    samples_in_units : array-like
        Array of parameter samples converted to the requested units.
    """
    param_key = get_parameter_key(param,p)
    
    if param in ['ecc', 'varpi', 'density']:
        par1, par2 = df[param_key[0]].values, df[param_key[1]].values
        
        if param=='ecc':
            samples = np.sqrt(par1**2 + par2**2)
        if param=='varpi':
            samples = np.arctan2(par1,par2) * (180 / np.pi)  
        if param=='density':
            samples = par1/(par2**3)

    else:
        samples = df[param_key].values

    mstar = df['mass_star_m_sun'].values
    rstar = df['radius_star_r_sun'].values

    samples_in_units = convert_units(samples,param, mstar, rstar, units)

    return samples_in_units
        

def get_samples(data, param, nb_planets, units='star'):
    """
    Extracts posterior samples for a given parameter for all planets
    from either a dictionary of dataframes or a single dataframe. 
    Converts the samples into the specified units.

    Parameters:
    -----------
    data : dict or pandas.DataFrame
        - If dictionary, each item contains the samples (for all planets) of a unique analysis.
        - If a single dataframe, contains the samples for all planets in a single analysis.
    param : str
        Parameter to extract.
    nb_planets : int
        Number of planets to analyse.
    units : str, optional (default is 'star')
        Unit system to convert the mass/radius/density samples ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    samples_dict : dict or numpy.ndarray
        - If input is a dictionary, returns a dictionary where each value is a 2D array of samples for each planet.
        - If input is a dataframe, returns a 2D array of samples for each planet.
    
    """
    if isinstance(data, dict):
        samples_dict = dict.fromkeys(data.keys())
        for key, df in data.items():
            samples_all_planets = []
            for p in range(nb_planets):
                try:
                    samples_in_units = extract_samples(df, param,  p, units)
                    samples_all_planets.append(samples_in_units)
                except KeyError:
                    break
            samples_dict[key] = np.array(samples_all_planets)
        return samples_dict
    
    elif isinstance(data, pd.DataFrame):
        samples_all_planets = []
        for p in range(nb_planets):
            try:
                samples_in_units = extract_samples(data, param,  p, units)
                samples_all_planets.append(samples_in_units)
            except KeyError:
                break
        return np.array(samples_all_planets)
    
    else:
        raise ValueError("Input data must be a dictionary or a pandas DataFrame.")

def get_offset(p, x_samples, y_samples, keys=None):
    """
    Computes the offset scales for the median plots based on the range of x and y parameters.
    """
    if keys is not None:
        x_all = np.hstack([x_samples[key][p] for key in keys])
        y_all = np.hstack([y_samples[key][p] for key in keys])

    else:
        x_all = np.hstack([x_samples[p]])
        y_all = np.hstack([y_samples[p]])

    x_offset_scale = 0.015 * (x_all.max() - x_all.min()) if x_all.size else 0
    y_offset_scale = 0.015 * (y_all.max() - y_all.min()) if y_all.size else 0

    return x_offset_scale, y_offset_scale

def get_best_fit(samples):
    t0_quantiles = np.percentile(samples, [16, 50, 84])
    median = t0_quantiles[1]
    lower = t0_quantiles[1] - t0_quantiles[0]
    upper = t0_quantiles[2] - t0_quantiles[1]
    return median, lower, upper

def plot_median_samples(axs, x_values, y_values, x_offset, y_offset, offset_scale, label):
    """
    Plots the median values for the given x and y samples with 1-sigma error bars. 
    The data points are offset to avoid overlapping between multiple analyses.

    Parameters:
    -----------
    axs : matplotlib.axes.Axes
        The axis on which to plot the error bars.
    x_values : array-like
        The x parameter values to plot.
    y_values : array-like
        The y parameter values to plot.
    x_offset : float
        The offset to adjust the x-values in the plot.
    y_offset : float
        The offset to adjust the y-values in the plot.
    offset_scale : float
        The scale factor of the offsets scaled to the range of the parameters.
    label : str
        The label of datapoints for the legend.

    Returns:
    --------
    None
        Function directly plots the data points and errorbars to the axes object.
    """
    x_best_fit, x_lower, x_upper = get_best_fit(x_values)
    y_best_fit, y_lower, y_upper = get_best_fit(y_values)
    
    axs.errorbar(
        x_best_fit + offset_scale * x_offset, 
        y_best_fit + offset_scale * y_offset, 
        xerr=[[x_lower], [x_upper]], 
        yerr=[[y_lower], [y_upper]], 
        fmt='o', label=label, capsize=4, 
        capthick=2, elinewidth=1.5, 
        markersize=6, markeredgewidth=1 
    )  

def plot_histograms(data, param, nb_planets=1, units='star'):
    """
    Plots histograms of the specified planetary parameter. 
    If input is a dictionary of dataframes, creates separate subplots for each planet.
    If input is a single dataframe, all planets are plotted on the same subplot.

    Parameters:
    -----------
    data : dict or pandas.DataFrame
        - If dictionary, each key corresponds to a different analysis, and each value is a dataframe of samples.
        - If a single dataframe, contains the posterior samples for multiple planets from a unique analysis.
    param : str
        The parameter to plot.
    nb_planets : int, optional (default is 1)
        The number of planets to analyse.
    units : str, optional (default is 'star')
        Unit system to convert the mass/radius/density samples ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    None
        Function creates and displays histograms for each planet.
    """
    if isinstance(data, dict):
        fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
        axs = np.atleast_1d(axes)
        x_samples = get_samples(data, param, nb_planets, units)

        for p in range(nb_planets):
            ax = axs[p]
            for key in data.keys():
                if p < len(x_samples[key]):
                    x = x_samples[key][p]
                    ax.hist(x, bins=50, alpha=0.5, label=key)
            ax.set_xlabel(get_labels(param, units))
            ax.legend()
            ax.set_title(f'planet {p}')

    elif isinstance(data, pd.DataFrame):
        fig, ax = plt.subplots(figsize=(6, 5))
        x_samples = get_samples(data, param, nb_planets, units)

        for p in range(nb_planets):
            x = x_samples[p]
            ax.hist(x, bins=50, alpha=0.5, label=f"planet {p+1}")
        ax.set_xlabel(get_labels(param, units))
        ax.legend()

    plt.tight_layout()
    plt.show()


def plot_samples(data, x_param, y_param, nb_planets=1, units='star', plot_median=False):
    """
    Plots the posterior samples for two planetary parameters for each planet.
    If input is a dictionary of dataframes, creates separate subplots for each planet.
    If input is a single dataframe, all planets are plotted on the same subplot.

    Parameters:
    -----------
    data : dict or pandas.DataFrame
        - If dictionary, each key corresponds to a different analysis, and each value is a dataframe of samples.
        - If a single dataframe, contains the posterior samples for multiple planets from a unique analysis.
    x_param : str
        The parameter to plot on the x-axis.
    y_param : str
        The parameter to plot on the y-axis.
    nb_planets : int, optional (default is 1)
        The number of planets to analyse.
    units : str, optional (default is 'star')
        Unit system to convert the mass/radius/density samples ('star', 'sun', 'earth', 'jup', or 'SI').
    plot_median : bool, optional (default is False)
        Whether to plot the median values and error bars instead of scatter points.

    Returns:
    --------
    None
        Function creates scatter plots for the given parameters for each planet.
    """
    if isinstance(data, dict):
        fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
        axs = np.atleast_1d(axes)

        x_samples = get_samples(data, x_param, nb_planets, units)
        y_samples = get_samples(data, y_param, nb_planets, units)

        for p in range(nb_planets):
            ax = axs[p]
            if plot_median:
                x_offset, y_offset = get_offset(p, x_samples, y_samples, data.keys())

            for i, key in enumerate(data.keys()):
                if p < len(x_samples[key]) and p < len(y_samples[key]):
                    x = x_samples[key][p]
                    y = y_samples[key][p]

                    if plot_median:
                        offset_scale = (i - len(data) / 2)
                        plot_median_samples(ax, x, y, x_offset, y_offset, offset_scale, key)
                    else:
                        ax.scatter(x, y, alpha=0.1, label=key)

            ax.set_xlabel(get_labels(x_param, units))
            ax.set_ylabel(get_labels(y_param, units))
            ax.legend()
            ax.set_title(f'planet {p}')

    elif isinstance(data, pd.DataFrame):
        fig, ax = plt.subplots(figsize=(6, 5))
        x_samples = get_samples(data, x_param, nb_planets, units)
        y_samples = get_samples(data, y_param, nb_planets, units)

        for p in range(nb_planets):
            x = x_samples[p]
            y = y_samples[p]
            label = f"planet {p+1}"

            if plot_median:
                x_offset, y_offset = np.median(x), np.median(y)
                plot_median_samples(ax, x, y, x_offset, y_offset, 0, label)
            else:
                ax.scatter(x, y, alpha=0.1, label=label)

        ax.set_xlabel(get_labels(x_param, units))
        ax.set_ylabel(get_labels(y_param, units))
        ax.legend()

    plt.tight_layout()
    plt.show()


def compare_period_ratios(df_dict, nb_planets=2):
    """
    Compares and plots the period ratios between consecutive planets.

    Parameters:
    -----------
    df_dict : dict
        Dictionary where each key corresponds to a unique analysis, and each value is a dataframe of samples.
    nb_planets : int, optional (default is 2)
        The number of planets to compare (must be at least 2).

    Returns:
    --------
    None
        Function creates histograms of period ratios between consecutive planets.
    """
    assert nb_planets >= 2, 'Need at least two planets to compare consecutively.'

    fig, axes = plt.subplots(nb_planets - 1, 1, figsize=(6, 5 * (nb_planets - 1)))
    axs = np.atleast_1d(axes)
    samples = get_samples(df_dict, 'period', nb_planets)

    for p in range(nb_planets - 1):
        ax = axs[p]
        for key in df_dict.keys():
            if p + 1 < len(samples[key]):
                val_p, val_p1 = samples[key][p], samples[key][p + 1]
                comparison = val_p1 / val_p
                label = f"P_{p+1} / P_{p}"
                ax.hist(comparison, bins=50, alpha=0.5, label=key)
        ax.set_xlabel(label + f" [{get_labels('period', units='star')}]")
        ax.legend()

    plt.tight_layout()
    plt.show()

def plot_consecutive_planets(df_dict, param, nb_planets=2, units='star'):
    """
    Plots the specified parameter for consecutive planets for each analysis.

    Parameters:
    -----------
    df_dict : dict
        Dictionary where each key corresponds to a unique analysis, and each value is a dataframe of samples.
    param : str
        The parameter to plot.
    nb_planets : int, optional (default is 2)
        The number of planets to compare (must be at least 2).
    units : str, optional (default is 'star')
        Unit system to convert the mass/radius/density samples ('star', 'sun', 'earth', 'jup', or 'SI').

    Returns:
    --------
    None
        Function creates scatter plots comparing consecutive planets for the given parameter.
    """
    assert nb_planets >= 2, 'Need at least two planets to compare consecutively.'

    fig, axes = plt.subplots(nb_planets - 1, 1, figsize=(6, 5 * (nb_planets - 1)))
    axs = np.atleast_1d(axes)
    samples = get_samples(df_dict, param, nb_planets, units)

    for p in range(nb_planets - 1):
        ax = axs[p]
        for i, key in enumerate(df_dict.keys()):
            if p + 1 < len(samples[key]):
                x, y = samples[key][p], samples[key][p + 1]
                ax.scatter(x, y, alpha=0.1, label=key)
        ax.set_xlabel(f'{get_labels(param, units)} (planet {p})')
        ax.set_ylabel(f'{get_labels(param, units)} (planet {p+1})')
        ax.set_title(f'{param}_{p+1} vs {param}_{p}')
        ax.legend()

    plt.tight_layout()
    plt.show()