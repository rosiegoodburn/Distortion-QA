"""
Crops a 3-D MRI volume to the phantom grid region.
Python equivalent of MATLAB's croppingTool.m.
"""

import numpy as np
from scipy.ndimage import uniform_filter1d
from skimage import morphology, restoration


def _smooth(x: np.ndarray, n: int) -> np.ndarray:
    """Moving-average smoothing (equivalent to MATLAB smooth())."""
    return uniform_filter1d(x.astype(float), size=n, mode='nearest')


def _gradient1d(x: np.ndarray) -> np.ndarray:
    """Central-difference gradient (like MATLAB gradient())."""
    return np.gradient(x.astype(float))


def _imflatfield(img: np.ndarray, sigma: float, mask: np.ndarray = None) -> np.ndarray:
    """
    Replicates MATLAB imflatfield: multiplicative flat-field correction.
    B = A * mean(A) / shading
    where shading is a Gaussian-smoothed version of A.
    """
    from scipy.ndimage import gaussian_filter
    A = img.astype(float)
    filter_size = 2 * int(np.ceil(2 * sigma)) + 1
    truncate = (filter_size - 1) / (2.0 * sigma)

    shading = gaussian_filter(A, sigma=sigma, truncate=truncate, mode='reflect')

    if mask is None:
        mean_val = np.nanmean(A)
    else:
        mean_val = np.nanmean(A[mask.astype(bool)])

    B = A * mean_val / (shading + 1e-12)
    B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)

    if mask is not None:
        B[~mask.astype(bool)] = A[~mask.astype(bool)]

    return B


def _bwareaopen(mask: np.ndarray, min_size: int) -> np.ndarray:
    return morphology.remove_small_objects(mask.astype(bool), min_size=min_size)


def cropping_tool(mri_old: np.ndarray) -> np.ndarray:
    """
    Crop a 3-D MRI volume to the phantom grid.

    Parameters
    ----------
    mri_old : 3-D float array (slices × rows × cols or similar)

    Returns
    -------
    mri_old_crop2 : 3-D cropped array with the thin dimension first
    """
    around_grid_crop_pix = 3
    top_btm_grid_crop_pix = 3

    # Ensure non-square dims are not first
    if mri_old.shape[1] == mri_old.shape[2]:
        mri_old = np.transpose(mri_old, (2, 1, 0))

    mri_old_sum1 = mri_old.sum(axis=0)   # sum over first dim → 2-D
    mri_old_sum2 = mri_old.sum(axis=1)   # sum over second dim → 2-D

    test1 = mri_old_sum1.max()
    test2 = mri_old_sum2.max()

    def _crop_2d_top(mri_old_sum):
        """Crop rows/cols from a top-view sum projection."""
        crop_loc_rows = _gradient1d(_smooth(mri_old_sum.sum(axis=0), 5))
        crop_rows_1 = int(np.argmax(crop_loc_rows[:-3])) + around_grid_crop_pix
        crop_rows_2 = int(np.argmin(crop_loc_rows[:-3])) - around_grid_crop_pix

        crop_loc_cols = _gradient1d(_smooth(mri_old_sum.sum(axis=1), 5))
        crop_cols_1 = int(np.argmax(crop_loc_cols)) + around_grid_crop_pix
        crop_cols_2 = int(np.argmin(crop_loc_cols)) - around_grid_crop_pix

        return crop_rows_1, crop_rows_2, crop_cols_1, crop_cols_2

    if test1 < test2:
        # ── Branch 1: top view from axis-0 sum ───────────────────────────────
        mri_old_sum = mri_old_sum1
        r1, r2, c1, c2 = _crop_2d_top(mri_old_sum)
        mri_old_crop = mri_old[:, c1:c2, r1:r2]

        # Side projection
        mri_old_sum2 = mri_old_crop.sum(axis=2)
        mask = (mri_old_sum2 >= 0.2 * mri_old_sum2.max()).astype(bool)
        if mask.sum() / mask.size > 0.3:
            mask = (mri_old_sum2 >= 0.5 * mri_old_sum2.max()).astype(bool)
        mask = _bwareaopen(mask, 2000)

        mri_old_sum2 = _imflatfield(mri_old_sum2, 15, mask)

        freq1 = _smooth(np.abs(_gradient1d(mri_old_sum2.sum(axis=0))), 50)
        freq2 = _smooth(np.abs(_gradient1d(mri_old_sum2.sum(axis=1))), 50)

        if freq1.max() > freq2.max():
            x_acc = np.stack([_smooth(np.abs(_gradient1d(mri_old_sum2[:, a])), 50)
                              for a in range(mri_old_sum2.shape[1])], axis=1)
            x_ = x_acc.sum(axis=1)
        else:
            x_acc = np.stack([_smooth(np.abs(_gradient1d(mri_old_sum2[a, :])), 50)
                              for a in range(mri_old_sum2.shape[0])], axis=0)
            x_ = x_acc.sum(axis=1)

        crop_locs = np.where(x_ > x_.max() / 2)[0]
        lo = int(crop_locs.min()) + top_btm_grid_crop_pix
        hi = int(crop_locs.max()) - top_btm_grid_crop_pix
        mri_old_crop2 = mri_old_crop[lo:hi, :, :]

    else:
        # ── Branch 2: top view from axis-1 sum ───────────────────────────────
        mri_old_sum = mri_old_sum2
        r1, r2, c1, c2 = _crop_2d_top(mri_old_sum)
        mri_old_crop = mri_old[c1:c2, :, r1:r2]

        # Side projection
        mri_old_sum2 = mri_old_crop.sum(axis=2)
        mask = (mri_old_sum2 >= 0.2 * mri_old_sum2.max()).astype(bool)
        if mask.sum() / mask.size > 0.3:
            mask = (mri_old_sum2 >= 0.5 * mri_old_sum2.max()).astype(bool)
        mask = _bwareaopen(mask, 2000)

        mri_old_sum2 = _imflatfield(mri_old_sum2, 15, mask)

        freq1 = _smooth(np.abs(_gradient1d(mri_old_sum2.sum(axis=0))), 50)
        freq2 = _smooth(np.abs(_gradient1d(mri_old_sum2.sum(axis=1))), 50)

        if freq1.max() > freq2.max():
            x_acc = np.stack([_smooth(np.abs(_gradient1d(mri_old_sum2[:, a])), 50)
                              for a in range(mri_old_sum2.shape[1])], axis=1)
            x_ = x_acc.sum(axis=0)
        else:
            x_acc = np.stack([_smooth(np.abs(_gradient1d(mri_old_sum2[a, :])), 50)
                              for a in range(mri_old_sum2.shape[0])], axis=0)
            x_ = x_acc.sum(axis=0)

        crop_locs = np.where(x_ > x_.max() / 2)[0]
        lo = int(crop_locs.min()) + top_btm_grid_crop_pix
        hi = int(crop_locs.max()) - top_btm_grid_crop_pix
        mri_old_crop2 = mri_old_crop[:, lo:hi, :]

    # ── Ensure the thin dimension is axis-0 ──────────────────────────────────
    dims = mri_old_crop2.shape
    min_dim = int(np.argmin(dims))
    if min_dim != 0:
        order = list(range(3))
        order[0], order[min_dim] = order[min_dim], order[0]
        mri_old_crop2 = np.transpose(mri_old_crop2, order)

    # ── Trim to at most 24 slices in first dimension ──────────────────────────
    if mri_old_crop2.shape[0] > 24:
        cut  = mri_old_crop2.shape[0] - 24
        cut1 = cut // 2
        cut2 = cut - cut1
        mri_old_crop2 = mri_old_crop2[cut1:mri_old_crop2.shape[0] - cut2, :, :]

    return mri_old_crop2