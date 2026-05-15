"""
Combines multiple 3-D distortion maps by taking the element-wise maximum
(ignoring NaN values where any map has data).
"""

import numpy as np


def combine_arrays(*arrays: np.ndarray) -> np.ndarray:
    """
    Parameters
    ----------
    *arrays : one or more np.ndarray of identical shape

    Returns
    -------
    np.ndarray
        Element-wise NaN-aware maximum across all input arrays.
    """
    if len(arrays) == 0:
        raise ValueError('At least one input required.')

    stacked = np.stack([a for a in arrays], axis=-1)  # (..., n)
    return np.nanmax(stacked, axis=-1)