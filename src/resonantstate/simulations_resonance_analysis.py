import numpy as np
import pandas as pd

def get_periods(data, row_number):
    if row_number == -1:
        row_idx = data.index[-1]
    else:
        row_idx = row_number

    # Get rows
    row_data = data[data.columns[data.columns.str.contains('period')]].iloc[row_idx].sort_values()
    row_index = np.vectorize(lambda x: int(x.split('_')[-1]))(np.array(row_data.index))
    row_values = np.array(row_data.values)
    
    # Filter nans
    nonan_inds = ~np.isnan(row_values)
    row_index = row_index[nonan_inds]
    row_values = row_values[nonan_inds]

    return row_values

def reorder_data(data, row_number):
    """
    Reorders a dataframe of posterior samples to be in increasing order of the periods.

    Args:
        data (pd.DataFrame): DataFrame containing posterior samples from one analysis.
        row_number (int): Index of the row based on which the columns will be reordered. If -1, the last row is used.

    Returns:
        pd.DataFrame: copy of the dataframe with the planet columns reordered in increasing period.
    """

    # Identify columns with planet data and their ids
    planet_cols = []
    planet_ids = set()
    param_names = [
        'mean_longitude_deg_',
        'period_days_',
        'k_',
        'h_',
        'inclination_deg_',
        'longitude_of_ascending_node_deg_',
        'mass_planet_star_ratio_',
        'radius_planet_star_ratio_'
    ]

    for col in data.columns:
        parts = col.split('_')
        if parts[-1].isdigit():
            planet_cols.append(col)
            planet_ids.add(int(parts[-1]))

    # Get the periods for each planet at the specified row
    period_cols = [f'period_days_{p_id}' for p_id in planet_ids if f'period_days_{p_id}' in data.columns]
    row_idx = data.index[-1] if row_number == -1 else row_number
    row_periods = data.loc[row_idx, period_cols]
    valid_periods = row_periods.dropna()

    # Sort the planet ids based on their periods at the specified row
    original_ids = [int(col.split('_')[-1]) for col in valid_periods.index]
    sorted_ids = [x for _, x in sorted(zip(valid_periods.values, original_ids))]

    # Map the old column names to new ones based on sorted ids
    column_map = {}
    for new_id, old_id in enumerate(sorted_ids):
        for param in param_names:
            old_col = param + str(old_id)
            new_col = param + str(new_id)
            if old_col in data.columns:
                column_map[old_col] = new_col

    # Separate planet columns and global columns and create new dataframe
    planet_cols_sorted = list(column_map.keys())
    reordered_df = data[planet_cols_sorted].rename(columns=column_map).copy()
    global_cols = [col for col in data.columns if col not in planet_cols_sorted]
    reordered_df[global_cols] = data[global_cols]
    if 'sample_index' in reordered_df.columns:
        reordered_df = reordered_df[['sample_index'] + [col for col in reordered_df.columns if col != 'sample_index']]


    return reordered_df




def get_nearest_resonance(period_ratio, second_order = False, kmax=12, difference_order = 0.5):
    """Simple function to get the nearest resonance for a pair of planets based on their period ratio.

    Args:
        period_ratio (float): period ratio between a pair of planets in the system.
        second_order (bool, optional): if second order resonances are allowed. Defaults to False.
        kmax (int, optional): maximum k value to be checked. Defaults to 12.
        difference_order (float, optional): penalty to be applied to width of second-order resonances. Defaults to 0.5.

    Returns:
        tuple: resonance k, order and Delta value associated to the closest resonance
    """
    
    # Get closest first order resonance
    ks = np.arange(1, kmax)
    possible_Delta = period_ratio * ks / (ks + 1) - 1
    k = ks[np.argmin(np.abs(possible_Delta))]
    min_Delta = np.min(np.abs(possible_Delta))
    
    # Get closest second order resonance
    if second_order:
        possible_Delta = period_ratio * ks / (ks + 2) - 1
        k2 = ks[np.argmin(np.abs(possible_Delta))]
        min_Delta2 = np.min(np.abs(possible_Delta))
        
        # If second order matches better first order, return it
        if min_Delta2 < difference_order * min_Delta:
            return k2, 2, min_Delta2
    
    return k, 1, min_Delta


def get_near_resonant_pairs(data, row_number, max_Delta = 1e-2):
    """Gets the near-resonant pairs for a given time step of the simulations

    Args:
        data (pd.DataFrame): time evolution of simulation of a planetary system
        row_number (int): index of the desired row to analyse
        max_Delta (float, optional): maximum value allowed for Delta to consider pair in resonance. Defaults to 1e-2.

    Returns:
        list[tuple]: list of resonance pairs with corresponding k value and order, respectively.
    """
    
    # Get rows
    row_idx = data.index[-1] if row_number == -1 else row_number
    row_data = data[data.columns[data.columns.str.contains('period')]].iloc[row_idx].sort_values()
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
    k_resonances, order, Deltas_resonances  = np.vectorize(lambda x: get_nearest_resonance(x))(period_ratios)
    order = np.array(order)
    k_resonances = np.array(k_resonances)
    
    # Filter resonances fulfilling the criterion to be considered in resonance
    inds = np.where(Deltas_resonances < max_Delta)
    resonances = list(zip(pairs_planets[inds], k_resonances[inds], order[inds]))
    
    return resonances
