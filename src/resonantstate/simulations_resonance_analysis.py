import numpy as np


def get_nearest_resonance(period_ratios, second_order = False, kmax=12, difference_order = 0.5):
    """Simple function to get the nearest resonance.

    Args:
        period_ratios (list): list of period ratios between the planets in the system.
        second_order (bool, optional): if second order resonances are allowed. Defaults to False.
        kmax (int, optional): maximum k value to be checked. Defaults to 12.

    Returns:
        tuple: resonance k, order and delta value associated to the closest resonance
    """
    
    # Get closest first order resonance
    ks = np.arange(2, kmax)
    possible_deltas = period_ratios * (ks - 1) / ks - 1
    k = ks[np.argmin(np.abs(possible_deltas))]
    min_delta = np.min(np.abs(possible_deltas))
    
    # Get closest second order resonance
    if second_order:
        possible_deltas = period_ratios * (ks - 2) / ks - 1
        k2 = ks[np.argmin(np.abs(possible_deltas))]
        min_delta2 = np.min(np.abs(possible_deltas))
        
        # If second order matches better first order, return it
        if min_delta2 < difference_order * min_delta:
            return k2, 2, min_delta2
    
    return k, 1, min_delta
