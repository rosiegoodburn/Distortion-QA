"""
Builds the ideal (undistorted) position map for all phantom grid cross-points.
"""

import numpy as np


def build_ideal_pos_map(orien_mm, grid_pos_mm):
    """
    Parameters
    ----------
    orien_mm     : list/array of 3 signed ints – phantom orientation vector
    grid_pos_mm  : list/array of 3 floats – grid centre position in mm

    Returns
    -------
    ideal_sideways : 21×21 float array (NaN where phantom rod absent)
    ideal_updown   : 21×21 float array (NaN where phantom rod absent)
    """
    orien_mm    = list(orien_mm)
    grid_pos_mm = np.asarray(grid_pos_mm, dtype=float)

    n = 21; c = 11; step = 14.25
    i_idx, j_idx = np.meshgrid(np.arange(1, n + 1), np.arange(1, n + 1),
                                indexing='ij')
    updown   = step * (i_idx - c)
    sideways = step * (j_idx - c)

    exclude = int(np.argmin(np.abs(orien_mm)))   # 0-based index

    dims = [0, 1, 2]
    dims.pop(exclude)

    # Orientation-specific sign correction
    if exclude == 1:    # horizontal
        pass
    elif exclude == 2:  # vertical-1
        updown = -updown
    else:               # vertical-2
        pass

    ideal_sideways = sideways + grid_pos_mm[dims[0]]
    ideal_updown   = updown   + grid_pos_mm[dims[1]]

    # ── Phantom hole mask ────────────────────────────────────────────────────
    # 1-based row/col pairs of missing rods
    mask = np.ones((21, 21), dtype=float)
    row = [7, 8, 7, 8, 9, 0, 1, 2, 8, 9, 10, 11, 18, 19, 20,
           0, 1, 2, 9, 10, 11, 18, 19, 20, 12, 13, 12, 13]
    col = [7, 7, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
           10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 12, 13, 13]
    for r, c_ in zip(row, col):
        mask[r, c_] = np.nan

    # ── 24 possible phantom orientations ─────────────────────────────────────
    base_mask = mask.copy()

    # horizontal (exclude Y = index 1)
    G = []; M = []

    mH = np.flipud(base_mask)
    G.append([-2, -1,  2]); M.append(mH.copy())
    G.append([ 2, -1,  2]); M.append(np.rot90(M[-1]))
    G.append([ 2, -1, -2]); M.append(np.rot90(M[-1]))
    G.append([-2, -1, -2]); M.append(np.rot90(M[-1]))

    G.append([ 2,  1,  2]); M.append(np.fliplr(mH))
    G.append([ 2,  1, -2]); M.append(np.rot90(M[-1]))
    G.append([-2,  1, -2]); M.append(np.rot90(M[-1]))
    G.append([-2,  1,  2]); M.append(np.rot90(M[-1]))

    # vertical-1 (exclude Z = index 2)
    mV1 = np.flipud(base_mask)
    G.append([-2, -2, -1]); M.append(mV1.copy())
    G.append([ 2, -2, -1]); M.append(np.rot90(M[-1]))
    G.append([ 2,  2, -1]); M.append(np.rot90(M[-1]))
    G.append([-2,  2, -1]); M.append(np.rot90(M[-1]))

    G.append([ 2, -2,  1]); M.append(np.fliplr(mV1))
    G.append([ 2,  2,  1]); M.append(np.rot90(M[-1]))
    G.append([-2,  2,  1]); M.append(np.rot90(M[-1]))
    G.append([-2, -2,  1]); M.append(np.rot90(M[-1]))

    # vertical-2 (exclude X = index 0)
    mV2 = np.rot90(base_mask)
    G.append([-1, -2,  2]); M.append(mV2.copy())
    G.append([-1,  2,  2]); M.append(np.rot90(M[-1]))
    G.append([-1,  2, -2]); M.append(np.rot90(M[-1]))
    G.append([-1, -2, -2]); M.append(np.rot90(M[-1]))

    G.append([ 1,  2,  2]); M.append(np.fliplr(mV2))
    G.append([ 1,  2, -2]); M.append(np.rot90(M[-1]))
    G.append([ 1, -2, -2]); M.append(np.rot90(M[-1]))
    G.append([ 1, -2,  2]); M.append(np.rot90(M[-1]))

    new_mask = None
    for g, m in zip(G, M):
        if g == orien_mm:
            new_mask = m
            break

    if new_mask is None:
        raise ValueError(f'Unrecognised orientation: {orien_mm}')

    ideal_sideways = ideal_sideways * new_mask
    ideal_updown   = ideal_updown   * new_mask

    return ideal_sideways, ideal_updown