"""
Shift utilities – equivalents of MATLAB's shift1px_zero_fill and shift2d_fill.
"""

import numpy as np


def shift1px_zero_fill(X: np.ndarray, d_row: int, d_col: int) -> np.ndarray:
    """
    Shift a 2-D boolean array by (d_row, d_col) with zero-fill at boundaries.

    Parameters
    ----------
    X     : 2-D bool array
    d_row : row shift (positive = down)
    d_col : column shift (positive = right)

    Returns
    -------
    Y : 2-D bool array, same shape as X
    """
    Y = np.zeros(X.shape, dtype=bool)

    r_src = np.arange(X.shape[0]) - d_row
    c_src = np.arange(X.shape[1]) - d_col

    r_ok = (r_src >= 0) & (r_src < X.shape[0])
    c_ok = (c_src >= 0) & (c_src < X.shape[1])

    Y[np.ix_(r_ok, c_ok)] = X[np.ix_(r_src[r_ok], c_src[c_ok])]
    return Y


def shift2d_fill(X: np.ndarray, d_row: int, d_col: int, fill_val) -> np.ndarray:
    """
    Shift a 2-D array by (d_row, d_col) with fill_val at boundaries.

    Parameters
    ----------
    X        : 2-D array (any numeric dtype)
    d_row    : row shift (positive = down)
    d_col    : column shift (positive = right)
    fill_val : scalar fill value (e.g. NaN or 0)

    Returns
    -------
    Y : 2-D array, same shape as X
    """
    Y = np.full(X.shape, fill_val, dtype=float)

    r_src = np.arange(X.shape[0]) - d_row
    c_src = np.arange(X.shape[1]) - d_col

    r_ok = (r_src >= 0) & (r_src < X.shape[0])
    c_ok = (c_src >= 0) & (c_src < X.shape[1])

    Y[np.ix_(r_ok, c_ok)] = X[np.ix_(r_src[r_ok], c_src[c_ok])]
    return Y