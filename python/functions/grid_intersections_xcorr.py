"""
Detect grid intersection points using cross-shaped template matching
(cross-correlation).
Python equivalent of MATLAB's gridIntersections_xcorr.m
"""

import numpy as np
from skimage import morphology, filters, measure
from skimage.morphology import reconstruction, local_maxima
from skimage.segmentation import clear_border
from scipy.ndimage import convolve, binary_erosion, binary_dilation
import cv2

def _imextendedmax(img: np.ndarray, h: float) -> np.ndarray:
    seed = img - h
    dilated = reconstruction(seed, img, method='dilation')
    # Regional maxima of the h-maxima transform: pixels strictly above reconstruction
    hmax = dilated  # this is imhmax output
    # imregionalmax: connected regions where hmax equals its regional max
    return local_maxima(hmax, connectivity=2)  # connectivity=2 → 8-connected in 2D


def grid_intersections_xcorr(img: np.ndarray,
                             L: int,
                             w: int,
                             angles=0) -> np.ndarray:
    """
    Detect phantom grid intersection (cross) points via cross-correlation.

    Parameters
    ----------
    img    : 2-D float or uint array (grayscale MRI slice)
    L      : half-length of each arm in pixels (~0.6-0.9 x grid spacing)
    w      : arm half-width in pixels (~grid line half-thickness)
    angles : scalar or iterable of rotation angles in degrees (default 0)

    Returns
    -------
    pts : Nx2 float array of [x, y] pixel coordinates (sub-pixel accuracy)
    """
    if np.isscalar(angles):
        angles = [angles]

    # Prepare image: invert so dark grid lines become bright
    I = img.astype(float)
    I = 1.0 - ((I - I.min()) / (I.max() - I.min() + 1e-12))

    Rmax = np.zeros(I.shape, dtype=float)

    sz = 2 * L + 1
    c  = L  # 0-based centre index

    for ang in angles:
        # Build cross kernel
        K = np.zeros((sz, sz), dtype=float)
        K[c, :]  = 1.0   # horizontal arm
        K[:, c]  = 1.0   # vertical arm

        if w > 0:
            selem = morphology.diamond(w)
            pad = max(0, w)
            K_padded = np.pad(K, pad)
            K_padded = binary_dilation(K_padded.astype(bool), selem).astype(float)
            K = K_padded[pad:pad+sz, pad:pad+sz]

        if ang != 0:
            M_rot = cv2.getRotationMatrix2D((c, c), ang, 1.0)
            K = cv2.warpAffine(K, M_rot, (sz, sz), flags=cv2.INTER_LINEAR)

        # Smooth, zero-mean, unit-norm normalise
        K = filters.gaussian(K, sigma=max(1, w))
        K = K - K.mean()
        norm_k = np.linalg.norm(K)
        if norm_k > 0:
            K /= norm_k

        R = convolve(I, K, mode='reflect')
        Rmax = np.maximum(Rmax, R)

    # Mask: suppress the bright rod regions themselves
    se3 = morphology.disk(3)
    se8 = morphology.disk(8)
    mask = I > filters.threshold_otsu(I)
    mask2 = binary_dilation(binary_erosion(mask, se3), se8)
    mask2 = ~mask2
    
    # Normalise masked response to [0, 1]  (equivalent to MATLAB mat2gray)
    Rmax_m = Rmax * mask2
    rng = Rmax_m.max() - Rmax_m.min()
    if rng > 0:
        Rmax_m = (Rmax_m - Rmax_m.min()) / rng

    # Peak detection: imextendedmax -> imclearborder -> bwareaopen -> imclose
    M_08 = _imextendedmax(Rmax_m, 0.08)
    M_08[0, :]  = 1; M_08[-1:, :] = 1; M_08[:, 0]  = 1; M_08[:, -1:] = 1

    M_02 = _imextendedmax(Rmax_m, 0.02)
    M_02[0, :]  = 1; M_02[-1:, :] = 1; M_02[:, 0]  = 1; M_02[:, -1:] = 1

    M_08 = clear_border(M_08); M_02 = clear_border(M_02)
    M_08 = morphology.remove_small_objects(M_08, min_size=1)
    M_02 = morphology.remove_small_objects(M_02, min_size=1)
    M_08 = morphology.binary_closing(M_08, morphology.disk(max(1, round(w / 2))))
    M_02 = morphology.binary_closing(M_02, morphology.disk(max(1, round(w / 2))))

    M = M_08 | M_02

    # Weighted centroids for sub-pixel accuracy
    labeled = measure.label(M)
    regions = measure.regionprops(labeled, intensity_image=Rmax_m)

    if not regions:
        return np.empty((0, 2))

    pts = np.array([r.weighted_centroid for r in regions])
    pts = pts[:, ::-1]   # (row, col) -> (x=col, y=row)
    return pts