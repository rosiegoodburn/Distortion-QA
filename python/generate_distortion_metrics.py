"""
Generates distortion metrics from measurements of the 2D UAL (Uniformity
and Linearity) phantom.

2026, Dr. Rosie Goodburn, Royal Marsden Hospital, UK

Phantom set-up:
  - Phantom must be completely full of water, including knob compartment
  - Ensure phantom is not twisted/oblique wrt XYZ
  - Appropriate setup should ensure that the isocentre lie within the
    phantom grid

Phantom can be rotated by 90°, flipped, or translated:
  - Horizontal, perpendicular-to-Y orientations for coronal 2D metrics
  - Vertical-1, perpendicular-to-Z orientations for axial 2D metrics
  - Vertical-2, perpendicular-to-X orientations for sagittal 2D metrics
"""

import os
from pathlib import Path
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

# Add functions directory to path
script_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_folder, 'functions'))

from generate_distortion_grid import generate_distortion_grid
from combine_arrays import combine_arrays
from calc_distortion_metrics import calc_distortion_metrics

# ── Measurements to include ───────────────────────────────────────────────────
data_root = Path('/Users/rgoodburn/Documents/MATLAB/Distortion')

# direc = data_root / 'scans 1'
# folders = ['3_fl3d_vibe_PE_RL_AXIAL_horizontal', '7_fl3d_vibe_PE_AP_SAGITTAL_vertical1', '10_fl3d_vibe_PE_HF_CORONAL_vertical2',]

# direc = data_root / 'scans 2'
# folders = ['4_fl3d_vibe_PE_RL_AXIAL_horizontal', '10_fl3d_vibe_PE_AP_SAGITTAL_vertical1', '13_fl3d_vibe_PE_HF_CORONAL_vertical2',]

# direc = data_root / 'scans 3'
# folders = ['11_fl3d_vibe_PE_RL_AXIAL_horizontal', '19_fl3d_vibe_PE_AP_SAGITTAL_vertical1', '22_fl3d_vibe_PE_HF_CORONAL_vertical2',]

# direc = data_root / 'scans 4'
# folders = ['3_fl3d_vibe_PE_RL_AXIAL_horizontal', '8_fl3d_vibe_PE_AP_SAGITTAL_vertical1', '13_fl3d_vibe_PE_HF_CORONAL_vertical2',]

direc = data_root / 'scans 5'
folders = ['3_fl3d_vibe_PE_RL_AXIAL_horizontal', '7_fl3d_vibe_PE_AP_SAGITTAL_vertical1', '10_fl3d_vibe_PE_HF_CORONAL_vertical2',]

scan_names = ['horizontal', 'vertical1', 'vertical2']

# ── Calculate 2D-distortion magnitude map ────────────────────────────────────
n = len(folders)
dis_map_mag_2D  = [None] * n
mri_grids       = [None] * n
pts_xs          = [None] * n
pts_ys          = [None] * n

for i in range(n):
    dis_map_mag_2D[i], mri_grids[i], pts_xs[i], pts_ys[i] = generate_distortion_grid(
        direc, folders[i], interp_phantom_grid_gaps=True)

# ── Save QA images ────────────────────────────────────────────────────────────
for i in range(n):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(mri_grids[i], cmap='gray')
    valid = (pts_xs[i] != 0) & (pts_ys[i] != 0)
    ax.plot(pts_xs[i][valid], pts_ys[i][valid], 'r+', markersize=10, markeredgewidth=1.5)
    ax.set_title(f'{scan_names[i].capitalize()} — detected grid points\n{folders[i]}')
    ax.axis('off')
    plt.tight_layout()
    img_path = os.path.join(direc, f'qa_grid_points_{scan_names[i]}.png')
    fig.savefig(img_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved QA image: {img_path}')

# ── Estimate 3D-distortion magnitude map ─────────────────────────────────────
# Assumption: Z² ≈ (X² + Y²) / 2  →  3D ≈ sqrt(2D² + (2D/√2)²)
dis_map_mag_3D_est = [None] * n
for i in range(n):
    A = dis_map_mag_2D[i]
    dis_map_mag_3D_est[i] = np.sqrt(A**2 + (A / np.sqrt(2))**2)

# ── Combine measurements (max where overlapping) ─────────────────────────────
combined_dis_map_mag_3D_est = combine_arrays(*dis_map_mag_3D_est)

# ── Derive 2D metrics ─────────────────────────────────────────────────────────
coords = np.arange(-250, 251)          # -250 … 250 inclusive (501 points)
X, Y, Z = np.meshgrid(coords, coords, coords, indexing='xy')

radius_50  = 50
radius_100 = 100

sphere_50  = (X**2 + Y**2 + Z**2) <= radius_50**2   # 100 mm DSV
sphere_100 = (X**2 + Y**2 + Z**2) <= radius_100**2  # 200 mm DSV

dis_map_mag_2D_100 = [None] * n
dis_map_mag_2D_200 = [None] * n
metrics_2D_100mm   = [None] * n
metrics_2D_200mm   = [None] * n

for i in range(n):
    m100 = sphere_50  * dis_map_mag_2D[i]
    m200 = sphere_100 * dis_map_mag_2D[i]

    m100 = m100.astype(float); m100[m100 == 0] = np.nan
    m200 = m200.astype(float); m200[m200 == 0] = np.nan

    dis_map_mag_2D_100[i] = m100
    dis_map_mag_2D_200[i] = m200

    metrics_2D_100mm[i] = calc_distortion_metrics(m100)
    metrics_2D_200mm[i] = calc_distortion_metrics(m200)

# ── Derive 3D metrics ─────────────────────────────────────────────────────────
m3d_100 = sphere_50  * combined_dis_map_mag_3D_est
m3d_200 = sphere_100 * combined_dis_map_mag_3D_est

m3d_100 = m3d_100.astype(float); m3d_100[m3d_100 == 0] = np.nan
m3d_200 = m3d_200.astype(float); m3d_200[m3d_200 == 0] = np.nan

metrics_3D_est_100mm = calc_distortion_metrics(m3d_100)
metrics_3D_est_200mm = calc_distortion_metrics(m3d_200)

# ── Save CSV ──────────────────────────────────────────────────────────────────
csv_path = os.path.join(direc, 'distortion_metrics.csv')
rows = []
for i, folder in enumerate(folders):
    rows.append({'scan': folder, 'volume': '2D_100mm', **metrics_2D_100mm[i]})
    rows.append({'scan': folder, 'volume': '2D_200mm', **metrics_2D_200mm[i]})
rows.append({'scan': 'combined', 'volume': '3D_est_100mm', **metrics_3D_est_100mm})
rows.append({'scan': 'combined', 'volume': '3D_est_200mm', **metrics_3D_est_200mm})

fieldnames = ['scan', 'volume', 'mean_distortion', 'std_distortion',
              'max_distortion', 'P95_distortion']
with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print(f'Results saved to {csv_path}')