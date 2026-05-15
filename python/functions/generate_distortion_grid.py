"""
Generates a 2-D distortion magnitude map from a folder of DICOM files.
Python equivalent of MATLAB's generateDistortionGrid.m

Rosie Goodburn, Royal Marsden Hospital
"""

import os
import numpy as np
import pydicom
from scipy.ndimage import gaussian_filter
from skimage.restoration import estimate_sigma

from find_phantom_pos import find_phantom_pos
from cropping_tool import cropping_tool
from grid_intersections_xcorr import grid_intersections_xcorr
from group_grid import group_grid_rows, group_grid_cols
from build_ideal_pos_map import build_ideal_pos_map
from transform_to_match import transform_to_match
from shift_utils import shift1px_zero_fill, shift2d_fill
from correct_pos_grids import correct_pos_grids
from add_distortion_map_to_new_grid import add_distortion_map_to_new_grid


def _natsort_key(s):
    """Natural sort key so DICOM filenames sort correctly (e.g. IM-0001-0010)."""
    import re
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r'(\d+)', s)]


def _imflatfield(img: np.ndarray, sigma: float) -> np.ndarray:
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
    mean_val = np.nanmean(A)
    B = A * mean_val / (shading + 1e-12)
    B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)
    return B


def _fillmissing_nearest(arr: np.ndarray) -> np.ndarray:
    """Fill NaN values using nearest non-NaN neighbour (1-D)."""
    arr = arr.copy()
    nans = np.isnan(arr)
    if not nans.any() or nans.all():
        return arr
    idx = np.arange(len(arr))
    arr[nans] = np.interp(idx[nans], idx[~nans], arr[~nans],
                          left=arr[~nans][0], right=arr[~nans][-1])
    return arr


def _remove_spurious_border_points(pts: np.ndarray, row_y: np.ndarray, 
                                    col_x: np.ndarray, min_spacing: float) -> np.ndarray:
    """
    Remove points belonging to the outermost row or column if it is closer
    than min_spacing pixels to the next one inward.
    """
    if len(pts) == 0:
        return pts

    keep = np.ones(len(pts), dtype=bool)
    row_y = np.array(sorted(row_y))
    col_x = np.array(sorted(col_x))

    # ── Rows ──────────────────────────────────────────────────────────────────
    if len(row_y) >= 2:
        if row_y[1] - row_y[0] < min_spacing:
            keep &= pts[:, 1] > row_y[0] + min_spacing / 2
        if row_y[-1] - row_y[-2] < min_spacing:
            keep &= pts[:, 1] < row_y[-1] - min_spacing / 2

    # ── Cols ──────────────────────────────────────────────────────────────────
    if len(col_x) >= 2:
        if col_x[1] - col_x[0] < min_spacing:
            keep &= pts[:, 0] > col_x[0] + min_spacing / 2
        if col_x[-1] - col_x[-2] < min_spacing:
            keep &= pts[:, 0] < col_x[-1] - min_spacing / 2

    return pts[keep]


def generate_distortion_grid(direc_mr: str,
                             folder_mr: str,
                             interp_phantom_grid_gaps: bool = True) -> np.ndarray:
    """
    Parameters
    ----------
    direc_mr                 : path to the parent scan directory
    folder_mr                : name of the scan subfolder
    interp_phantom_grid_gaps : whether to interpolate over missing grid rods

    Returns
    -------
    new_grid_mag : (501, 501, 501) float array – distortion magnitude in mm
    """

    # ── Load DICOM volume ─────────────────────────────────────────────────────
    dicom_folder = os.path.join(direc_mr, folder_mr, 'DICOM')
    dicom_files  = sorted(
        [f for f in os.listdir(dicom_folder) if f.endswith('.dcm')],
        key=_natsort_key)

    info = pydicom.dcmread(os.path.join(dicom_folder, dicom_files[0]))
    vol  = np.zeros((info.Rows, info.Columns, len(dicom_files)), dtype=np.uint16)
    for i, fname in enumerate(dicom_files):
        ds = pydicom.dcmread(os.path.join(dicom_folder, fname))
        vol[:, :, i] = ds.pixel_array

    # Build a plain dict of the metadata needed downstream
    info_mr = {
        'InPlanePhaseEncodingDirection': info.InPlanePhaseEncodingDirection,
        'ImageOrientationPatient':       np.asarray(info.ImageOrientationPatient, float),
        'ImagePositionPatient':          np.asarray(info.ImagePositionPatient, float),
        'PixelSpacing':                  np.asarray(info.PixelSpacing, float),
        'SliceThickness':                float(info.SliceThickness),
    }

    # ── Find phantom orientation and grid centre ───────────────────────────────
    orien_mm, grid_pos_mm = find_phantom_pos(vol, info_mr)

    pix_sz = info_mr['PixelSpacing'][0]
    slc_sz = info_mr['SliceThickness']

    # ── Crop volume to phantom grid ───────────────────────────────────────────
    mri_cropped = cropping_tool(vol.astype(float))

    # ── Bias-field correction and projection ──────────────────────────────────
    mri_cropped_br = np.zeros_like(mri_cropped, dtype=float)
    for i in range(mri_cropped.shape[0]):
        mri_cropped_br[i] = _imflatfield(mri_cropped[i], 15)
    mri_grid = mri_cropped_br.sum(axis=0)   # 2-D projection

    # ── First pass: find intersection points ──────────────────────────────────
    pts = grid_intersections_xcorr(mri_grid, L=8, w=1, angles=0)
    _, row_idx, row_y = group_grid_rows(pts, tol=1)
    _, col_idx, col_x = group_grid_cols(pts, tol=1)

    pts = _remove_spurious_border_points(pts, row_y, col_x, min_spacing=7.0)

    # redo grouping on cleaned pts
    _, row_idx, row_y = group_grid_rows(pts, tol=1)
    _, col_idx, col_x = group_grid_cols(pts, tol=1)

    pts_x = pts[:, 0].copy()
    pts_y = pts[:, 1].copy()
    pts_x[row_idx == 0] = 0
    pts_y[row_idx == 0] = 0
    pts_x[col_idx == 0] = 0
    pts_y[col_idx == 0] = 0

    pts_ = np.column_stack([pts_x, pts_y])
    pts_ = pts_[~np.any(pts_ == 0, axis=1)]
    
    pts_x = pts_[:, 0].copy()
    pts_y = pts_[:, 1].copy()

    # Convert pixel positions to mm
    if mri_grid.shape[0] > mri_grid.shape[1]:
        pts_[:, 0] *= slc_sz
        pts_[:, 1] *= pix_sz
    else:
        pts_[:, 0] *= pix_sz
        pts_[:, 1] *= slc_sz

    # ── Build coarse 21×21 grid ───────────────────────────────────────────────
    step = 14.25
    pts_grid = pts_ / step
    start = np.array([pts_grid[:, 0].min(), pts_grid[:, 1].min()])
    pts_grid = np.round(pts_grid + (1 - start)).astype(int)

    grid_m = np.full((21, 21), np.nan)
    for i in range(len(pts_grid)):
        r, c = pts_grid[i, 1] - 1, pts_grid[i, 0] - 1   # 0-based
        if 0 <= r < 21 and 0 <= c < 21:
            grid_m[r, c] = 1.0

    # ── Ideal position map ────────────────────────────────────────────────────
    ideal_sideways, ideal_updown = build_ideal_pos_map(orien_mm, grid_pos_mm)

    # ── Detect missing row/column ─────────────────────────────────────────────
    no_row_1 = np.all(np.isnan(grid_m[0,  :]))
    no_row_2 = np.all(np.isnan(grid_m[-1, :]))
    no_row_3 = np.all(np.isnan(grid_m[:,  0]))
    no_row_4 = np.all(np.isnan(grid_m[:, -1]))
    row_missing = no_row_1 or no_row_2 or no_row_3 or no_row_4

    side_missing = None
    if row_missing:
        test = (~np.isnan(ideal_sideways)).astype(float)
        if test[:, :9].sum() == 179:
            side_missing = 1
            ideal_sideways[:, 0] = np.nan; ideal_updown[:, 0] = np.nan
        elif test[:9, :].sum() == 179:
            side_missing = 2
            ideal_sideways[0, :] = np.nan; ideal_updown[0, :] = np.nan
        elif test[:, 12:].sum() == 179:
            side_missing = 3
            ideal_sideways[:, -1] = np.nan; ideal_updown[:, -1] = np.nan
        elif test[12:, :].sum() == 179:
            side_missing = 4
            ideal_sideways[-1, :] = np.nan; ideal_updown[-1, :] = np.nan

    # ── Transform image to match ideal orientation ────────────────────────────
    mri_grid_new = transform_to_match(mri_grid, grid_m, ideal_sideways)

    # ── Second pass: intersection detection on aligned image ──────────────────
    pts = grid_intersections_xcorr(mri_grid_new, L=8, w=1, angles=0)
    _, row_idx, row_y = group_grid_rows(pts, tol=1)
    _, col_idx, col_x = group_grid_cols(pts, tol=1)

    pts = _remove_spurious_border_points(pts, row_y, col_x, min_spacing=7.0)

    # redo grouping on cleaned pts
    _, row_idx, row_y = group_grid_rows(pts, tol=1)
    _, col_idx, col_x = group_grid_cols(pts, tol=1)

    pts_x = pts[:, 0].copy()
    pts_y = pts[:, 1].copy()
    pts_x[row_idx == 0] = 0
    pts_y[row_idx == 0] = 0
    pts_x[col_idx == 0] = 0
    pts_y[col_idx == 0] = 0

    pts_ = np.column_stack([pts_x, pts_y])
    pts_ = pts_[~np.any(pts_ == 0, axis=1)]

    pts_x = pts_[:, 0].copy()
    pts_y = pts_[:, 1].copy()

    if mri_grid_new.shape[0] > mri_grid_new.shape[1]:
        pts_[:, 0] *= slc_sz
        pts_[:, 1] *= pix_sz
    else:
        pts_[:, 0] *= pix_sz
        pts_[:, 1] *= slc_sz

    # ── Build measured position grids ─────────────────────────────────────────
    pts_grid = pts_ / step
    start = np.array([pts_grid[:, 0].min(), pts_grid[:, 1].min()])
    pts_grid = np.round(pts_grid + (1 - start)).astype(int)

    grid1  = np.zeros((21, 21))
    grid2  = np.zeros((21, 21))
    grid_m = np.full((21, 21), np.nan)

    for i in range(len(pts_grid)):
        r, c = pts_grid[i, 1] - 1, pts_grid[i, 0] - 1
        if 0 <= r < 21 and 0 <= c < 21:
            grid_m[r, c] = 1.0
            grid1[r, c]  = pts_[i, 0]
            grid2[r, c]  = pts_[i, 1]

    grid1 *= grid_m
    grid2 *= grid_m

    # ── Shift grids if a row is missing ──────────────────────────────────────
    if row_missing:
        test_ideal = (~np.isnan(ideal_sideways)).astype(int)
        test_grid  = (~np.isnan(grid1)).astype(int)

        shifts = [(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)]
        best_score = np.inf
        best_shift = (0, 0)
        for s in shifts:
            shifted = shift1px_zero_fill(test_grid, s[0], s[1])
            score   = np.count_nonzero(np.logical_xor(shifted, test_ideal))
            if score < best_score:
                best_score = score
                best_shift = s

        d_row, d_col = best_shift
        grid1 = shift2d_fill(grid1, d_row, d_col, np.nan)
        grid2 = shift2d_fill(grid2, d_row, d_col, np.nan)

    # ── Interpolate over phantom-grid gaps ────────────────────────────────────
    if interp_phantom_grid_gaps:
        for i in range(21):
            grid1[:, i]          = _fillmissing_nearest(grid1[:, i])
            ideal_sideways[:, i] = _fillmissing_nearest(ideal_sideways[:, i])
            grid2[i, :]          = _fillmissing_nearest(grid2[i, :])
            ideal_updown[i, :]   = _fillmissing_nearest(ideal_updown[i, :])

        if row_missing and side_missing is not None:
            if side_missing == 1:
                grid1[:, 0] = np.nan; grid2[:, 0] = np.nan
                ideal_sideways[:, 0] = np.nan; ideal_updown[:, 0] = np.nan
            elif side_missing == 2:
                grid1[0, :] = np.nan; grid2[0, :] = np.nan
                ideal_sideways[0, :] = np.nan; ideal_updown[0, :] = np.nan
            elif side_missing == 3:
                grid1[:, -1] = np.nan; grid2[:, -1] = np.nan
                ideal_sideways[:, -1] = np.nan; ideal_updown[:, -1] = np.nan
            elif side_missing == 4:
                grid1[-1, :] = np.nan; grid2[-1, :] = np.nan
                ideal_sideways[-1, :] = np.nan; ideal_updown[-1, :] = np.nan

    # ── Correct grid positioning ──────────────────────────────────────────────
    grid1, grid2 = correct_pos_grids(
        orien_mm, grid1, grid2, ideal_sideways, ideal_updown, row_missing)

    # ── Distortion maps ───────────────────────────────────────────────────────
    distortion_map1 = ideal_sideways - grid1
    distortion_map2 = ideal_updown   - grid2

    # ── Resample onto 501×501×501 volume ─────────────────────────────────────
    new_grid = add_distortion_map_to_new_grid(
        grid_pos_mm, orien_mm, distortion_map1, distortion_map2)

    new_grid = np.nan_to_num(new_grid, nan=0.0)

    # Magnitude: RSS across the 3 scanner dimensions (axis 0)
    new_grid_mag = np.sqrt(np.sum(new_grid**2, axis=0))

    return new_grid_mag, mri_grid_new, pts_x, pts_y