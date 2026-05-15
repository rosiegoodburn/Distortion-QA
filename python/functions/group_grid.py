"""
Group detected intersection points into rows and columns using
single-linkage hierarchical clustering.
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster


def _cluster_1d(values: np.ndarray, tol: float):
    """Return integer cluster labels for 1-D coordinate array."""
    if len(values) == 1:
        return np.array([1])
    Z = linkage(values[:, None], method='single')
    return fcluster(Z, t=tol, criterion='distance')


def group_grid_rows(pts: np.ndarray, tol: float = 1.0, min_count: int = 8):
    """
    Group points into rows based on their y-coordinate.

    Parameters
    ----------
    pts       : N×2 array [x, y]
    tol       : max y-distance within a row (pixels)
    min_count : discard rows with fewer points

    Returns
    -------
    groups : list of index arrays (each sorted left→right)
    row_idx : N-length label array (0 = discarded)
    row_y   : median y for each surviving row (sorted top→bottom)
    """
    y = pts[:, 1].astype(float)
    g = _cluster_1d(y, tol)

    unique_g = np.unique(g)
    groups = []
    row_y  = []
    for k in unique_g:
        idx = np.where(g == k)[0]
        order = np.argsort(pts[idx, 0])   # left→right
        groups.append(idx[order])
        row_y.append(np.median(y[idx]))

    # Sort top→bottom
    ord_r  = np.argsort(row_y)
    groups = [groups[i] for i in ord_r]
    row_y  = [row_y[i]  for i in ord_r]

    # Filter by size
    keep   = [i for i, gr in enumerate(groups) if len(gr) >= min_count]
    groups = [groups[i] for i in keep]
    row_y  = np.array([row_y[i] for i in keep])

    # Map to labels
    row_idx = np.zeros(len(y), dtype=int)
    for k, gr in enumerate(groups, start=1):
        row_idx[gr] = k

    return groups, row_idx, row_y


def group_grid_cols(pts: np.ndarray, tol: float = 1.0, min_count: int = 8):
    """
    Group points into columns based on their x-coordinate.

    Parameters
    ----------
    pts       : N×2 array [x, y]
    tol       : max x-distance within a column (pixels)
    min_count : discard columns with fewer points

    Returns
    -------
    groups : list of index arrays (each sorted top→bottom)
    col_idx : N-length label array (0 = discarded)
    col_x   : median x for each surviving column (sorted left→right)
    """
    x = pts[:, 0].astype(float)
    g = _cluster_1d(x, tol)

    unique_g = np.unique(g)
    groups = []
    col_x  = []
    for k in unique_g:
        idx = np.where(g == k)[0]
        order = np.argsort(pts[idx, 1])   # top→bottom
        groups.append(idx[order])
        col_x.append(np.median(x[idx]))

    # Sort left→right
    ord_c  = np.argsort(col_x)
    groups = [groups[i] for i in ord_c]
    col_x  = [col_x[i]  for i in ord_c]

    # Filter by size
    keep   = [i for i, gr in enumerate(groups) if len(gr) >= min_count]
    groups = [groups[i] for i in keep]
    col_x  = np.array([col_x[i] for i in keep])

    # Map to labels
    col_idx = np.zeros(len(x), dtype=int)
    for k, gr in enumerate(groups, start=1):
        col_idx[gr] = k

    return groups, col_idx, col_x
