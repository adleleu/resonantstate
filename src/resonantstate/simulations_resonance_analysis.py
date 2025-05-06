import numpy as np


def get_nearest_resonance(period_ratio, second_order = False, kmax=12, difference_order = 0.5):
    """Simple function to get the nearest resonance for a pair of planets based on their period ratio.

    Args:
        period_ratio (float): period ratio between a pair of planets in the system.
        second_order (bool, optional): if second order resonances are allowed. Defaults to False.
        kmax (int, optional): maximum k value to be checked. Defaults to 12.

    Returns:
        tuple: resonance k, order and delta value associated to the closest resonance
    """
    
    # Get closest first order resonance
    ks = np.arange(2, kmax)
    possible_delta = period_ratio * (ks - 1) / ks - 1
    k = ks[np.argmin(np.abs(possible_delta))]
    min_delta = np.min(np.abs(possible_delta))
    
    # Get closest second order resonance
    if second_order:
        possible_delta = period_ratio * (ks - 2) / ks - 1
        k2 = ks[np.argmin(np.abs(possible_delta))]
        min_delta2 = np.min(np.abs(possible_delta))
        
        # If second order matches better first order, return it
        if min_delta2 < difference_order * min_delta:
            return k2, 2, min_delta2
    
    return k, 1, min_delta


def get_resonance_pairs(data, row_number, max_delta = 3e-3):
    """Gets the resonance pairs for a given time step of the simulations

    Args:
        data (pd.DataFrame): time evolution of simulation of a planetary system
        row_number (int): index of the desired row to analyse
        max_delta (float, optional): maximum value allowed for delta to consider pair in resonance. Defaults to 3e-3.

    Returns:
        list[tuple]: list of resonance pairs with corresponding k value and order, respectively.
    """
    
    # Get rows
    row_data = data[data.columns[data.columns.str.contains('period')]].iloc[row_number].sort_values()
    row_index = np.vectorize(lambda x: int(x.split('_')[-1]))(np.array(row_data.index))
    row_values = np.array(row_data.values)
    
    # Filter nans
    nonan_inds = ~np.isnan(row_values)
    row_index = row_index[nonan_inds]
    row_values = row_values[nonan_inds]
    
    # Calculate Period Ratios and planet pairs
    period_ratios = row_values[1:] / row_values[:-1]
    pairs_planets = np.array(list(zip(row_index[:-1], row_index[1:])))
    
    # Compute closest resonances
    k_resonances, order, deltas_resonances  = np.vectorize(lambda x: get_nearest_resonance(x))(period_ratios)
    order = np.array(order)
    k_resonances = np.array(k_resonances)
    
    # Filter resonances fulfilling the criterion to be considered in resonance
    inds = np.where(deltas_resonances < max_delta)
    resonances = list(zip(pairs_planets[inds], k_resonances[inds], order[inds]))
    
    return resonances