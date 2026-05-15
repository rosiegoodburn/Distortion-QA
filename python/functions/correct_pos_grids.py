"""
Corrects the positioning of the measured grids so that the central region
matches the ideal (isocentre-aligned) grid positions.
"""

import numpy as np


def correct_pos_grids(orien_mm, grid1: np.ndarray, grid2: np.ndarray,
                      ideal_sideways: np.ndarray, ideal_updown: np.ndarray,
                      row_missing: bool = False):
    """
    Parameters
    ----------
    orien_mm       : list/array of 3 signed ints – phantom orientation vector
    grid1          : 2-D measured sideways positions
    grid2          : 2-D measured up-down positions
    ideal_sideways : 2-D ideal sideways positions
    ideal_updown   : 2-D ideal up-down positions
    row_missing    : bool (unused, kept for API compatibility)

    Returns
    -------
    grid1, grid2 : corrected position arrays
    """
    orien_mm = np.asarray(orien_mm)
    test = int(np.argmin(np.abs(orien_mm))) + 1  # 1-based like MATLAB

    if test == 3:   # vertical-1
        grid2 = -grid2
        grid2 = grid2 + -1 * (np.nanmin(grid2) + np.nanmax(grid2))

    # Find the isocentre row/col in the ideal grids
    min_sideways = np.nanmin(np.abs(ideal_sideways))
    _, min_sideways_col = np.where(np.abs(ideal_sideways) == min_sideways)

    min_updown = np.nanmin(np.abs(ideal_updown))
    min_updown_row, _ = np.where(np.abs(ideal_updown) == min_updown)

    row_c = int(min_updown_row[0])
    col_c = int(min_sideways_col[0])

    rows = slice(max(0, row_c - 3), min(ideal_sideways.shape[0], row_c + 4))
    cols = slice(max(0, col_c - 3), min(ideal_sideways.shape[1], col_c + 4))

    ideal_sw_mean  = np.nanmean(ideal_sideways[rows, cols])
    ideal_ud_mean  = np.nanmean(ideal_updown[rows, cols])
    grid1_mean     = np.nanmean(grid1[rows, cols])
    grid2_mean     = np.nanmean(grid2[rows, cols])

    grid1 = grid1 - (grid1_mean - ideal_sw_mean)
    grid2 = grid2 - (grid2_mean - ideal_ud_mean)

    return grid1, grid2