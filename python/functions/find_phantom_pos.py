"""
Determines the phantom orientation (orien_mm) and grid centre position
(grid_pos_mm) in scanner frame-of-reference from the MRI volume and DICOM
metadata.
"""

import numpy as np
from scipy.ndimage import binary_closing, binary_opening, binary_fill_holes
from skimage import morphology, filters, measure


def _binarize(volume: np.ndarray, threshold: float) -> np.ndarray:
    return volume > threshold


def _bwareaopen(mask: np.ndarray, min_size: int) -> np.ndarray:
    return morphology.remove_small_objects(mask.astype(bool), min_size=min_size)


def _sphere_selem(r: int):
    """3-D spherical structuring element of radius r."""
    d = 2 * r + 1
    c = r
    z, y, x = np.ogrid[:d, :d, :d]
    return ((x - c)**2 + (y - c)**2 + (z - c)**2) <= r**2


def find_phantom_pos(volume_mr: np.ndarray, info_mr: dict):
    """
    Parameters
    ----------
    volume_mr : 3-D uint16 array [rows, cols, slices]
    info_mr   : dict with DICOM header fields:
                  InPlanePhaseEncodingDirection, ImageOrientationPatient,
                  ImagePositionPatient, PixelSpacing, SliceThickness

    Returns
    -------
    orien_mm    : list of 3 signed ints (e.g. [-2, -1, 2])
    grid_pos_mm : 1-D float array of length 3
    """
    phase_direc = info_mr['InPlanePhaseEncodingDirection']
    if phase_direc == 'ROW':
        freq_direc = 1
    else:
        freq_direc = 2

    se2 = _sphere_selem(2)
    se4 = _sphere_selem(4)

    vol_f = volume_mr.astype(float)
    im_binary = binary_closing(
        binary_opening(_binarize(vol_f, 20), se2), se4)

    if im_binary.sum() > im_binary.size / 2:
        im_binary = binary_closing(
            binary_opening(_binarize(vol_f, 220), se2), se4)

    im_binary = _bwareaopen(im_binary, 50000)

    if freq_direc == 1:
        # top-down projection: sum over axis-1 (columns)
        from_top = np.stack(
            [im_binary[:, i, :].sum(axis=0) for i in range(im_binary.shape[1])],
            axis=0)  # shape: (cols, slices)

        knob_from_top = from_top.copy()
        knob_from_top[knob_from_top < 0.85 * knob_from_top.max()] = 0
        knob_bin = binary_closing(knob_from_top > 0, morphology.disk(2).reshape(1, -1, 1)[:, :, 0].astype(bool))
        # Simplified: just threshold + largest connected component
        from skimage.morphology import remove_small_objects, label
        knob_lab = label(knob_from_top > 0)
        if knob_lab.max() > 0:
            sizes = np.bincount(knob_lab.ravel())
            sizes[0] = 0
            knob_from_top = (knob_lab == sizes.argmax())
        rs, cs = np.where(knob_from_top)
        rowcolCOM1_knob = np.array([rs.mean(), cs.mean()])

        main_from_top = (from_top > 0)
        rs2, cs2 = np.where(main_from_top)
        rowcolCOM1_main = np.array([rs2.mean(), cs2.mean()])

        crop = int(round(rowcolCOM1_knob[0]))
        col_range = slice(max(0, crop - 40), min(im_binary.shape[1], crop + 41))
        from_side = np.stack(
            [im_binary[i, col_range, :].sum(axis=0) for i in range(im_binary.shape[0])],
            axis=0)

        # knob from side: suppress bright regions first, then binarize, open, area filter
        knob_fs = from_side.copy().astype(float)
        knob_fs[knob_fs > 0.6 * knob_fs.max()] = 0
        knob_fs_bin = knob_fs > filters.threshold_otsu(knob_fs) if knob_fs.max() > 0 else knob_fs > 0
        knob_fs_bin = binary_opening(knob_fs_bin, morphology.disk(4))
        knob_fs_bin = _bwareaopen(knob_fs_bin, 300)
        rs3, cs3 = np.where(knob_fs_bin)
        rowcolCOM2_knob = np.array([rs3.mean() if len(rs3) else 0,
                                    cs3.mean() if len(cs3) else 0])

        # main body from side
        main_fs = from_side.copy().astype(float)
        main_fs[main_fs < 0.85 * main_fs.max()] = 0
        main_fs_bin = binary_fill_holes(main_fs > 0)
        rs4, cs4 = np.where(main_fs_bin)
        rowcolCOM2_main = np.array([rs4.mean() if len(rs4) else 0,
                                    cs4.mean() if len(cs4) else 0])

        com_knob_pix = np.array([rowcolCOM2_knob[0], rowcolCOM1_knob[0],
                                  np.mean([rowcolCOM1_knob[1], rowcolCOM2_knob[1]])])
        com_main_pix = np.array([rowcolCOM2_main[0], rowcolCOM1_main[0],
                                  np.mean([rowcolCOM1_main[1], rowcolCOM2_main[1]])])
    else:
        from_top = np.stack(
            [im_binary[i, :, :].sum(axis=0) for i in range(im_binary.shape[0])],
            axis=0)

        knob_from_top = from_top.copy().astype(float)
        knob_from_top[knob_from_top < 0.85 * knob_from_top.max()] = 0
        from skimage.morphology import label
        knob_lab = label(knob_from_top > 0)
        if knob_lab.max() > 0:
            sizes = np.bincount(knob_lab.ravel())
            sizes[0] = 0
            knob_from_top = (knob_lab == sizes.argmax())
        rs, cs = np.where(knob_from_top)
        rowcolCOM1_knob = np.array([rs.mean(), cs.mean()])

        main_from_top = from_top > 0
        rs2, cs2 = np.where(main_from_top)
        rowcolCOM1_main = np.array([rs2.mean(), cs2.mean()])

        crop = int(round(rowcolCOM1_knob[0]))
        row_range = slice(max(0, crop - 40), min(im_binary.shape[0], crop + 41))
        from_side = np.stack(
            [im_binary[row_range, i, :].sum(axis=0) for i in range(im_binary.shape[1])],
            axis=0)

        knob_fs = from_side.copy().astype(float)
        knob_fs[knob_fs > 0.6 * knob_fs.max()] = 0
        knob_fs_bin = knob_fs > filters.threshold_otsu(knob_fs) if knob_fs.max() > 0 else knob_fs > 0
        knob_fs_bin = binary_opening(knob_fs_bin, morphology.disk(4))
        knob_fs_bin = _bwareaopen(knob_fs_bin, 300)
        rs3, cs3 = np.where(knob_fs_bin)
        rowcolCOM2_knob = np.array([rs3.mean() if len(rs3) else 0,
                                    cs3.mean() if len(cs3) else 0])

        main_fs = from_side.copy().astype(float)
        main_fs[main_fs < 0.85 * main_fs.max()] = 0
        main_fs_bin = binary_fill_holes(main_fs > 0)
        rs4, cs4 = np.where(main_fs_bin)
        rowcolCOM2_main = np.array([rs4.mean() if len(rs4) else 0,
                                    cs4.mean() if len(cs4) else 0])

        com_knob_pix = np.array([rowcolCOM1_knob[0], rowcolCOM2_knob[0],
                                  np.mean([rowcolCOM1_knob[1], rowcolCOM2_knob[1]])])
        com_main_pix = np.array([rowcolCOM1_main[0], rowcolCOM2_main[0],
                                  np.mean([rowcolCOM1_main[1], rowcolCOM2_main[1]])])

    # ── Convert pixel positions to mm ─────────────────────────────────────────
    orien_dcm = np.round(info_mr['ImageOrientationPatient']).astype(int)
    im_par_pos = np.asarray(info_mr['ImagePositionPatient'], dtype=float)
    pix_sz = float(info_mr['PixelSpacing'][0])
    slc_sz = float(info_mr['SliceThickness'])
    vox_sz = np.array([pix_sz, pix_sz, slc_sz])

    O2knob_mm = (com_knob_pix - vox_sz / 2) * vox_sz
    O2main_mm = (com_main_pix - vox_sz / 2) * vox_sz

    Opos_mm = np.zeros(3)
    if np.array_equal(orien_dcm, [1, 0, 0, 0, 1, 0]):
        Opos_mm[:] = im_par_pos
        knob_mm = _calc_pos(Opos_mm, O2knob_mm)
        main_mm = _calc_pos(Opos_mm, O2main_mm)
        knob_mm = np.array([knob_mm[1], knob_mm[0], knob_mm[2]])
        main_mm = np.array([main_mm[1], main_mm[0], main_mm[2]])
    elif np.array_equal(orien_dcm, [0, 1, 0, 0, 0, -1]):
        Opos_mm[:] = [im_par_pos[2], im_par_pos[1], im_par_pos[0]]
        knob_mm = _calc_pos(Opos_mm, O2knob_mm)
        main_mm = _calc_pos(Opos_mm, O2main_mm)
        knob_mm = np.array([knob_mm[2], knob_mm[1], knob_mm[0]])
        main_mm = np.array([main_mm[2], main_mm[1], main_mm[0]])
    elif np.array_equal(orien_dcm, [1, 0, 0, 0, 0, -1]):
        Opos_mm[:] = [im_par_pos[2], im_par_pos[0], im_par_pos[1]]
        knob_mm = _calc_pos(Opos_mm, O2knob_mm)
        main_mm = _calc_pos(Opos_mm, O2main_mm)
        knob_mm = np.array([knob_mm[1], knob_mm[2], knob_mm[0]])
        main_mm = np.array([main_mm[1], main_mm[2], main_mm[0]])
    else:
        raise ValueError('Unexpected orientation vector')

    # ── Determine orien_mm and grid_pos_mm ───────────────────────────────────
    main2knob_mm = knob_mm - main_mm
    min_i = int(np.argmin(np.abs(main2knob_mm)))

    # For near-zero components np.sign() returns 0 — default to +1 in that case
    orien = np.zeros(3, dtype=int)
    for k in range(3):
        s = int(np.sign(main2knob_mm[k]))
        orien[k] = s if s != 0 else 1

    orien_mm = (2 * orien).tolist()
    min_sign = int(np.sign(main2knob_mm[min_i]))
    orien_mm[min_i] = min_sign if min_sign != 0 else 1

    if min_i == 1:  # horizontal
        g2k = np.array([96.5037, 100.0875, 96.0920]) * orien
    elif min_i == 2:  # vertical-1
        g2k = np.array([96.3468,  98.9717,  99.5187]) * orien
    else:            # vertical-2
        g2k = np.array([100.2873,  99.2393,  95.4823]) * orien

    grid_pos_mm = knob_mm - g2k
    return orien_mm, grid_pos_mm


def _calc_pos(Opos_mm, O2x_mm):
    pos = np.zeros(3)
    for i in range(3):
        if Opos_mm[i] >= 0:
            pos[i] = Opos_mm[i] - O2x_mm[i]
        else:
            pos[i] = Opos_mm[i] + O2x_mm[i]
    return pos