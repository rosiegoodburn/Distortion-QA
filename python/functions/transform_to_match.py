"""
Finds the rotation/flip transformation that best aligns the detected grid
with the ideal grid, then applies it to the MRI image.
"""

import numpy as np


def transform_to_match(mri_grid: np.ndarray,
                       grid: np.ndarray,
                       ideal_sideways: np.ndarray) -> np.ndarray:
    """
    Parameters
    ----------
    mri_grid       : 2-D float array – summed MRI image of the phantom grid
    grid           : 2-D float array – binary (NaN/1) detected grid map
    ideal_sideways : 2-D float array – ideal sideways positions (NaN where absent)

    Returns
    -------
    mri_grid_new : 2-D float array – mri_grid after the best-matching transform
    """
    # Build 8 candidate orientations of the grid
    g1 = grid
    g2 = np.rot90(g1)
    g3 = np.rot90(g2)
    g4 = np.rot90(g3)
    g5 = np.fliplr(g1)
    g6 = np.rot90(g5)
    g7 = np.rot90(g6)
    g8 = np.rot90(g7)

    candidates = [g1, g2, g3, g4, g5, g6, g7, g8]

    # Mask of ideal positions (1 where not NaN)
    test = np.where(~np.isnan(ideal_sideways), 1.0, np.nan)

    # Score each orientation: nansum of (candidate + test)
    scores = [np.nansum(c + test) for c in candidates]

    best = int(np.argmax(scores))  # 0-indexed

    transforms = [
        lambda img: img,
        lambda img: np.rot90(img),
        lambda img: np.rot90(np.rot90(img)),
        lambda img: np.rot90(np.rot90(np.rot90(img))),
        lambda img: np.fliplr(img),
        lambda img: np.rot90(np.fliplr(img)),
        lambda img: np.rot90(np.rot90(np.fliplr(img))),
        lambda img: np.rot90(np.rot90(np.rot90(np.fliplr(img)))),
    ]

    return transforms[best](mri_grid)