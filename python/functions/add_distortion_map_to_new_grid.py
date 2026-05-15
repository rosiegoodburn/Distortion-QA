"""
Resamples the 2-D distortion maps onto a common 501×501×501 mm³ volume in
scanner frame-of-reference.
Python equivalent of MATLAB's addDistortionMapToNewGrid.m.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def add_distortion_map_to_new_grid(grid_pos_mm,
                                   orien_mm,
                                   distortion_map1: np.ndarray,
                                   distortion_map2: np.ndarray) -> np.ndarray:
    """
    Parameters
    ----------
    grid_pos_mm      : (3,) array – grid centre in scanner mm
    orien_mm         : list of 3 signed ints – phantom orientation
    distortion_map1  : 2-D array – sideways distortion (mm)
    distortion_map2  : 2-D array – up-down distortion (mm)

    Returns
    -------
    new_grid : (3, 501, 501, 501) float array
               Axis-0 = scanner dimension (X=0, Y=1, Z=2)
    """
    grid_pos_mm = np.asarray(grid_pos_mm, dtype=float)
    orien_mm    = np.asarray(orien_mm, dtype=float)

    coords = np.arange(-250, 251, dtype=float)            # 501 points
    Xq, Yq, Zq = np.meshgrid(coords, coords, coords, indexing='xy')

    new_grid = np.zeros((3, 501, 501, 501), dtype=float)

    test = int(np.argmin(np.abs(orien_mm))) + 1   # 1-based, as in MATLAB

    r, c = distortion_map1.shape
    # Pad distortion maps by 1 on each edge → (r+2, c+2)
    def _padded_vol(dm, slot, n_slices=3):
        """Return a (3, n_slices, r+2, c+2) volume with dm in the interior."""
        vol = np.full((3, n_slices, r + 2, c + 2), np.nan)
        vol[slot, :, 1:-1, 1:-1] = dm
        return vol

    step = 14.25

    if test == 2:   # ── horizontal orientation (freq enc = Y) ─────────────────
        dm_vol = np.full((3, 3, r + 2, c + 2), np.nan)
        dm_vol[0, :, 1:-1, 1:-1] = distortion_map1
        dm_vol[2, :, 1:-1, 1:-1] = distortion_map2
        dm_vol = np.transpose(dm_vol, (0, 1, 3, 2))   # permute [1 2 4 3]

        V1 = dm_vol[0]   # (3, c+2, r+2)
        V3 = dm_vol[2]

        X_pts = np.arange(grid_pos_mm[0] - step*11,
                          grid_pos_mm[0] + step*11 + step/2, step)
        Y_pts = np.arange(grid_pos_mm[1] - step,
                          grid_pos_mm[1] + step + step/2, step)
        Z_pts = np.arange(grid_pos_mm[2] - step*11,
                          grid_pos_mm[2] + step*11 + step/2, step)

        V1q = _interp_nearest(Y_pts, X_pts, Z_pts, V1, Xq, Yq, Zq)
        V3q = _interp_nearest(Y_pts, X_pts, Z_pts, V3, Xq, Yq, Zq)

        new_grid[0] = V1q
        new_grid[2] = V3q

    elif test == 3:   # ── vertical-1 orientation (freq enc = Z) ──────────────
        dm_vol = np.full((3, 3, r + 2, c + 2), np.nan)
        dm_vol[0, :, 1:-1, 1:-1] = distortion_map1
        dm_vol[1, :, 1:-1, 1:-1] = distortion_map2
        dm_vol = np.transpose(dm_vol, (0, 2, 3, 1))   # permute [1 3 4 2]

        V1 = dm_vol[0]
        V2 = dm_vol[1]

        X_pts = np.arange(grid_pos_mm[0] - step*11,
                          grid_pos_mm[0] + step*11 + step/2, step)
        Y_pts = np.arange(grid_pos_mm[1] + step*11,
                          grid_pos_mm[1] - step*11 - step/2, -step)
        Z_pts = np.arange(grid_pos_mm[2] + step,
                          grid_pos_mm[2] - step - step/2, -step)

        V1q = _interp_nearest(Y_pts, X_pts, Z_pts, V1, Xq, Yq, Zq)
        V2q = _interp_nearest(Y_pts, X_pts, Z_pts, V2, Xq, Yq, Zq)

        new_grid[0] = V1q
        new_grid[1] = V2q

    elif test == 1:   # ── vertical-2 orientation (freq enc = X) ──────────────
        dm_vol = np.full((3, 3, r + 2, c + 2), np.nan)
        dm_vol[1, :, 1:-1, 1:-1] = distortion_map1
        dm_vol[2, :, 1:-1, 1:-1] = distortion_map2
        dm_vol = np.transpose(dm_vol, (0, 3, 1, 2))   # permute [1 4 2 3]

        V2 = dm_vol[1]
        V3 = dm_vol[2]

        X_pts = np.arange(grid_pos_mm[0] - step,
                          grid_pos_mm[0] + step + step/2, step)
        Y_pts = np.arange(grid_pos_mm[1] - step*11,
                          grid_pos_mm[1] + step*11 + step/2, step)
        Z_pts = np.arange(grid_pos_mm[2] - step*11,
                          grid_pos_mm[2] + step*11 + step/2, step)

        V2q = _interp_nearest(Y_pts, X_pts, Z_pts, V2, Xq, Yq, Zq)
        V3q = _interp_nearest(Y_pts, X_pts, Z_pts, V3, Xq, Yq, Zq)

        new_grid[1] = V2q
        new_grid[2] = V3q

    return new_grid


def _interp_nearest(y_pts, x_pts, z_pts, V, Xq, Yq, Zq):
    """
    Nearest-neighbour 3-D interpolation.
    V has shape (ny, nx, nz) matching (y_pts, x_pts, z_pts).
    Returns an array of shape Xq.shape.
    """
    # Sort axes (scipy requires ascending coords)
    y_asc  = np.sort(y_pts);   iy = np.argsort(y_pts)
    x_asc  = np.sort(x_pts);   ix = np.argsort(x_pts)
    z_asc  = np.sort(z_pts);   iz = np.argsort(z_pts)

    V_sorted = V[np.ix_(iy, ix, iz)]

    interp = RegularGridInterpolator(
        (y_asc, x_asc, z_asc), V_sorted,
        method='nearest', bounds_error=False, fill_value=np.nan)

    pts = np.column_stack([Yq.ravel(), Xq.ravel(), Zq.ravel()])
    return interp(pts).reshape(Xq.shape)