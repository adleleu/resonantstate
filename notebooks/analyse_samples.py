import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from resonantstate.ell2SFM import *
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
    numplt.ndarray
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



def get_all_planets(data):
    planets_masterlist = []
    for df_dict in data:
        planets = df_dict["planets_list"]
        for planet in planets:
            if planet.lower() not in planets_masterlist:
                planets_masterlist.append(planet.lower())
    return planets_masterlist
        

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
    if type(dict_list) is not list:
        dict_list = [dict_list]

    all_planets = get_all_planets(dict_list)
    nb_planets = len(all_planets)
    
    fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
    axs = np.atleast_1d(axes)
    axis_dict = {all_planets[i]: axs[i] for i in range(nb_planets)}

    # Loop through each analysis and plot the histograms
    for df_dict in dict_list:
        analysis_id = df_dict['samples_name'].split('_')[-1]
        planets = df_dict['planets_list']
        df = df_dict['samples']

        for planet in planets:
            p_id = planets.index(planet)
            ax = axis_dict[planet.lower()]
            x = get_samples(df, param,  p_id, units)
            ax.hist(x, bins=50, alpha=0.5, label=f'analysis {analysis_id}')
            ax.set_xlabel(get_labels(param, units))
            ax.legend()
            ax.set_title(planet)

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
    if type(dict_list) is not list:
        dict_list = [dict_list]
    
    all_planets = get_all_planets(dict_list)
    nb_planets = len(all_planets)
    
    fig, axes = plt.subplots(nb_planets, 1, figsize=(6, 5 * nb_planets))
    axs = np.atleast_1d(axes)
    axis_dict = {all_planets[i]: axs[i] for i in range(nb_planets)}

    # Loop through each analysis and plot the samples
    for df_dict in dict_list:
        analysis_id = df_dict['samples_name'].split('_')[-1]
        df = df_dict['samples']
        planets = df_dict['planets_list']

        for planet in planets:
            p_id = planets.index(planet)
            ax = axis_dict[planet.lower()]
            ax.set_title(planet) 
            ax.set_xlabel(get_labels(x_param, units))
            ax.set_ylabel(get_labels(y_param, units))
            #if planet not in axis_dict:
            #    continue
            
            x = get_samples(df, x_param,  p_id, units)
            y = get_samples(df, y_param,  p_id, units)

            if (x.size == 0) or (y.size == 0):
                    continue
            ax.scatter(x, y, alpha=0.1, label=f'analysis {analysis_id}')
            ax.legend()   

    plt.tight_layout()
    plt.show()


def compare_period_ratios(dict_list, planet_pair):
    """
    Plots histograms of period ratios (P_{i+1} / P_i) between adjacent planets.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    planet_pair : list 
        List of planet pairs to be considered.
    """
    planet_pair_norm = [p.lower() for p in planet_pair]
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    for df_dict in dict_list:
        analysis_id = df_dict['samples_name'].split('_')[-1]
        df = df_dict['samples']
        planets = df_dict['planets_list']
        planets_norm = [p.lower() for p in planets]

        try:
            p_id0 = planets_norm.index(planet_pair_norm[0])
            p_id1 = planets_norm.index(planet_pair_norm[1])
        except ValueError:
            continue

        val_p0 = get_samples(df, 'period',  p_id0, units='star')
        val_p1 = get_samples(df, 'period',  p_id1, units='star')
        
        ax.hist(val_p1 / val_p0, bins=50, alpha=0.5, label=f'analysis {analysis_id}')
        ax.set_xlabel(f'Period ratio (P_{planet_pair[1]}/P_{planet_pair[0]})')
        ax.set_title(f'{planet_pair[1]} against {planet_pair[0]}')
        ax.legend()
    plt.tight_layout()
    plt.show()

def plot_adjacent_planets(dict_list, param, planet_pair, units='star'):
    """
    Plots scatter plots comparing a parameter between consecutive planets.

    Parameters:
    -----------
    dict_list : list of dict
        Each dictionary should contain 'sample_name', 'sample', and 'planets_list'.
    param : str
        Parameter to compare.
    planet_pair : list 
        List of planet pairs to be considered.
    units : str, optional (default is 'star')
        Target unit system ('star', 'sun', 'earth', 'jup', or 'SI'). 
        Applies only to parameters 'mass', 'radius', or 'density'.
    """
    planet_pair_norm = [p.lower() for p in planet_pair]
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    for df_dict in dict_list:
        analysis_id = df_dict['samples_name'].split('_')[-1]
        df = df_dict['samples']
        planets = df_dict['planets_list']
        planets_norm = [p.lower() for p in planets]

        try:
            p_id0 = planets_norm.index(planet_pair_norm[0])
            p_id1 = planets_norm.index(planet_pair_norm[1])
        except ValueError:
            continue

        x = get_samples(df, param,  p_id0, units)
        y = get_samples(df, param,  p_id1, units)  

        ax.scatter(x, y, alpha=0.1, label=f'analysis {analysis_id}')

        ax.set_xlabel(f'{get_labels(param, units)} ({planet_pair[0]})')
        ax.set_ylabel(f'{get_labels(param, units)} ({planet_pair[1]})')  
        ax.set_title(f'{planet_pair[1]} against {planet_pair[0]}')    
        ax.legend()   
    plt.tight_layout()
    plt.show()


def plot_ell2SFM_comparison(data_list, planet_pair, resonance, delta_lim=None, X_lim=None):
      """
      Compares the samples of different analyses in the phase space of the second fundamental model (SFM) of resonance.

      Parameters
      ----------
      data_list : list of dict
            List of dictionaries, each containing 'samples' (DataFrame) and 'analysis_id' (str).
      planet_pair : list or np.ndarray
            Pair of planets to be considered in the sample.
      resonance : int
            Resonance of the corresponding pair (p such that resonance is p:p+1).
      delta_lim : tuple, optional
            Lower and upper limits of the x-axis. If None, scales automatically with data.
      X_lim : tuple, optional
            Lower and upper limits of the y-axis. If None, scales automatically with data.
            
      Returns
      -------
      fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
            The figure and axes objects of the plot.
      """
      fig, ax = plt.subplots(figsize=(9, 9))

      cmap = plt.get_cmap('tab10')
      colors = [cmap(i % cmap.N) for i in range(len(data_list))]

      # Loop over each analysis in the list of dictionaries
      for i, data in enumerate(data_list):
            samples = data['samples']
            label = f'analysis {data["analysis_id"]}'
            print('Analysis ID:', data['analysis_id'])
            color = colors[i]

            if isinstance(planet_pair, np.ndarray):
                  pair = planet_pair.tolist()
            else:
                  pair = planet_pair

            plot_samples_SFM(fig, ax, samples, pair, resonance, colors=color, color_lim=None, label=label)
      plot_topology(ax, delta_lim, X_lim)

      plt.legend()
      plt.tight_layout()
      plt.show()
      return fig, ax




def plot_ell2SFM(data, colors=None, color_lim=None, delta_lim=None, X_lim=None):
      """
      Finds and plots the samples of near resonant planet pairs into the phase space of the second fundamental model (SFM) of resonance.

      Parameters
      ----------
      data : pandas.DataFrame or np.ndarray or dict
            Input data containing the posterior samples.
            - If DataFrame or np.ndarray: used directly as sample input.
            - If dict: must contain keys 'sample' (DataFrame) and 'samples_name' (str).
      colors : list of np.ndarray or list of str, optional 
            Numpy array used for colormapping the samples.
      color_lim : tuple, optional 
            Lower and upper limits for the color scale. If None, uses the min and max of the provided colormap array.
      delta_lim : tuple, optional
            Lower and upper limits of the x-axis. If None, scales automatically with data.
      X_lim : tuple, optional
            Lower and upper limits of the y-axis. If None, scales automatically with data.
      
      Returns
      -------
      fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
            The figure and axes objects of the plot.
      """

      if isinstance(data, dict):
            samples = data['samples']
            analysis_id = data['samples_name']
      elif isinstance(data, pd.DataFrame) or isinstance(data, np.ndarray):
            samples = data
      else:
            raise TypeError('Unsupported data type. Input has to be a pandas DataFrame, a numpy array, or a dictionary containing the "samples" key.')

      row = -1

      # Check if samples are ordered by increasing period
      if isinstance(samples, pd.DataFrame):
            periods_values = get_periods(samples, row)
            if not np.all(np.diff(periods_values) > 0):
                  samples = reorder_data(samples, row)

      # Get first-order resonant pairs
      #pairs_resonances = get_near_resonant_pairs(samples, row)
      #first_order_pairs = [(pair.tolist(), res, order) for pair, res, order in pairs_resonances if order == 1]

      #if not first_order_pairs:
      #      print('No first order pairs found.')
      #      return None, None

      # Create figure to plot samples if first order pairs are found
      fig, ax = py.subplots(1, 1, figsize=(9,9))
      if isinstance(data, dict):
            fig.suptitle(f'Analysis {analysis_id}', fontsize=16)

      # Print the pairs and their resonances
      #planet_pairs, p_indexes, orders = map(list, zip(*first_order_pairs))
      
      # Get first-order resonant pairs
      planet_pairs, p_indexes = get_p_by_pair(samples)
      
      print('Found', len(planet_pairs), 'first order pairs.')
      print('Pairs:', planet_pairs)
      print('Resonances:', p_indexes)

      # Use a colormap if colors is not provided
      if colors is None:
            cmap = py.get_cmap('tab10')
            colors = [cmap(i % cmap.N) for i in range(len(planet_pairs))]

      if any(isinstance(c, (np.ndarray, list)) for c in colors):
            flat_colors = np.concatenate([np.asarray(c).flatten() for c in colors])
            color_min = flat_colors.min() if color_lim is None else color_lim[0]
            color_max = flat_colors.max() if color_lim is None else color_lim[1]
            color_lim = (color_min, color_max)
            cbar=fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.hsv, 
                                                    norm=mpl.colors.Normalize(color_lim[0], color_lim[1])), 
                                                    ax=ax, aspect=40, pad=0.01)
            cbar.ax.tick_params()

      # Plot the samples
      plot_samples_SFM(fig, ax, samples, planet_pairs, p_indexes, colors, color_lim=color_lim, label='')
      plot_topology(ax, delta_lim, X_lim)
      py.legend()
      py.tight_layout()
      py.show()
      return fig, ax

