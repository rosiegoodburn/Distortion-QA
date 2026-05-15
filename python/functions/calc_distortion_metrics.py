"""
Calculates distortion metrics from a distortion magnitude map.
"""

import numpy as np


def calc_distortion_metrics(dis_map_mag: np.ndarray) -> dict:
    """
    Parameters
    ----------
    dis_map_mag : np.ndarray
        3-D distortion magnitude array (NaN where excluded).

    Returns
    -------
    dict with keys:
        mean_distortion, std_distortion, max_distortion, P95_distortion
    """
    flat = dis_map_mag.ravel()
    return {
        'mean_distortion': float(np.nanmean(flat)),
        'std_distortion':  float(np.nanstd(flat, ddof=0)),   # MATLAB std uses N-1; use ddof=1 to match exactly
        'max_distortion':  float(np.nanmax(flat)),
        'P95_distortion':  float(np.nanpercentile(flat, 95)),
    }