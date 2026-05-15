"""
XNAT runner for UAL distortion metrics.

Downloads 3 DICOM scan folders from XNAT, runs the distortion analysis,
and uploads the resulting CSV and QA images back to XNAT as session resources.

Usage:
    python xnat_runner.py \
        --url https://your-xnat.org \
        --user myusername \
        --password mypassword \
        --project PROJECT_ID \
        --subject SUBJECT_ID \
        --session SESSION_ID \
        --scan-axial    "3 fl3d_vibe_PE_RL_AXIAL_horizontal" \
        --scan-sagittal "8 fl3d_vibe_PE_AP_SAGITTAL_vertical1" \
        --scan-coronal  "13 fl3d_vibe_PE_HF_CORONAL_vertical2"
"""

import argparse
import csv
import os
import sys
import tempfile
import requests
import zipfile
import io

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — no display needed in Docker
import matplotlib.pyplot as plt
import numpy as np

# Add functions directory to path
script_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_folder, 'functions'))

from generate_distortion_grid import generate_distortion_grid
from combine_arrays import combine_arrays
from calc_distortion_metrics import calc_distortion_metrics


def parse_args():
    p = argparse.ArgumentParser(description='UAL phantom distortion analysis')
    p.add_argument('--url',           required=True,  help='XNAT base URL, e.g. https://your-xnat.org')
    p.add_argument('--user',          required=True,  help='XNAT username')
    p.add_argument('--password',      required=True,  help='XNAT password')
    p.add_argument('--project',       required=True,  help='XNAT project ID')
    p.add_argument('--subject',       required=True,  help='XNAT subject ID')
    p.add_argument('--session',       required=True,  help='XNAT session ID (experiment label)')
    p.add_argument('--scan-axial',    required=True,  help='Scan name for axial (horizontal) acquisition')
    p.add_argument('--scan-sagittal', required=True,  help='Scan name for sagittal (vertical-1) acquisition')
    p.add_argument('--scan-coronal',  required=True,  help='Scan name for coronal (vertical-2) acquisition')
    p.add_argument('--resource-name', default='UAL_DISTORTION_METRICS',
                   help='XNAT resource name for output files (default: UAL_DISTORTION_METRICS)')
    return p.parse_args()


def get_session(url: str, user: str, password: str) -> requests.Session:
    """Create an authenticated requests session."""
    s = requests.Session()
    s.auth = (user, password)
    r = s.get(f'{url}/data/auth', timeout=30)
    if r.status_code != 200:
        raise RuntimeError(f'XNAT authentication failed (HTTP {r.status_code}). Check your credentials.')
    print('Authenticated successfully.')
    return s


def get_scan_id(session: requests.Session, url: str,
                project: str, subject: str, experiment: str,
                scan_label: str) -> str:
    """Resolve a scan label to its XNAT scan ID."""
    r = session.get(
        f'{url}/data/projects/{project}/subjects/{subject}'
        f'/experiments/{experiment}/scans',
        params={'format': 'json'},
        timeout=30)
    r.raise_for_status()
    scans = r.json()['ResultSet']['Result']
    for scan in scans:
        if scan.get('series_description') == scan_label or scan.get('ID') == scan_label:
            return scan['ID']
    # Fall back to partial match
    for scan in scans:
        if scan_label in scan.get('series_description', '') or scan_label in scan.get('ID', ''):
            return scan['ID']
    available = [f"ID={s['ID']}, desc={s.get('series_description', '')}" for s in scans]
    raise ValueError(f'Scan "{scan_label}" not found in session {experiment}.\n'
                     f'Available scans:\n' + '\n'.join(available))


def download_scan(session: requests.Session, url: str,
                  project: str, subject: str, experiment: str,
                  scan_id: str, dest_folder: str):
    """Download all DICOM files for a scan into dest_folder/DICOM/."""
    dicom_dir = os.path.join(dest_folder, 'DICOM')
    os.makedirs(dicom_dir, exist_ok=True)

    r = session.get(
        f'{url}/data/projects/{project}/subjects/{subject}'
        f'/experiments/{experiment}/scans/{scan_id}/resources/DICOM/files',
        params={'format': 'zip'},
        timeout=300,
        stream=True)
    r.raise_for_status()

    # Extract only .dcm files, flattening any subfolder structure
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        for member in z.namelist():
            if member.lower().endswith('.dcm'):
                filename = os.path.basename(member)
                with z.open(member) as src, open(os.path.join(dicom_dir, filename), 'wb') as dst:
                    dst.write(src.read())

    n_files = len(os.listdir(dicom_dir))
    print(f'  Downloaded {n_files} DICOM files to {dicom_dir}')


def upload_file(session: requests.Session, url: str,
                project: str, subject: str, experiment: str,
                resource_name: str, file_path: str, content_type: str):
    """Upload a single file as a session resource on XNAT."""
    filename = os.path.basename(file_path)

    # Create resource if it doesn't exist
    session.put(
        f'{url}/data/projects/{project}/subjects/{subject}'
        f'/experiments/{experiment}/resources/{resource_name}',
        timeout=30)

    with open(file_path, 'rb') as f:
        r = session.put(
            f'{url}/data/projects/{project}/subjects/{subject}'
            f'/experiments/{experiment}/resources/{resource_name}/files/{filename}',
            data=f,
            headers={'Content-Type': content_type},
            timeout=60)
    r.raise_for_status()
    print(f'  Uploaded {filename} to XNAT resource {resource_name}')


def save_qa_image(mri_grid: np.ndarray, pts_x: np.ndarray, pts_y: np.ndarray,
                  title: str, output_path: str):
    """Save a QA image showing the grid with detected cross-points overlaid."""
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(mri_grid, cmap='gray')
    # Only plot valid (non-zero) points
    valid = (pts_x != 0) & (pts_y != 0)
    ax.plot(pts_x[valid], pts_y[valid], 'r+', markersize=10, markeredgewidth=1.5)
    ax.set_title(title, fontsize=12)
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved QA image: {output_path}')


def run_analysis(direc: str, folders: list, tmpdir: str) -> dict:
    """Run the distortion analysis and return all metrics plus QA image paths."""
    n = len(folders)
    scan_names = ['axial', 'sagittal', 'coronal']

    print('Generating 2D distortion maps...')
    dis_map_mag_2D = [None] * n
    qa_image_paths = []

    for i in range(n):
        print(f'  Processing scan {i+1}/{n}: {folders[i]}')
        mag, mri_grid_new, pts_x, pts_y = generate_distortion_grid(
            direc, folders[i], interp_phantom_grid_gaps=True)
        dis_map_mag_2D[i] = mag

        # Save QA image
        img_path = os.path.join(tmpdir, f'qa_grid_points_{scan_names[i]}.png')
        save_qa_image(mri_grid_new, pts_x, pts_y,
                      title=f'{scan_names[i].capitalize()} — detected grid points\n{folders[i]}',
                      output_path=img_path)
        qa_image_paths.append(img_path)

    print('Estimating 3D distortion maps...')
    dis_map_mag_3D_est = [None] * n
    for i in range(n):
        A = dis_map_mag_2D[i]
        dis_map_mag_3D_est[i] = np.sqrt(A**2 + (A / np.sqrt(2))**2)

    combined = combine_arrays(*dis_map_mag_3D_est)

    print('Calculating metrics...')
    coords = np.arange(-250, 251)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='xy')
    sphere_50  = (X**2 + Y**2 + Z**2) <= 50**2
    sphere_100 = (X**2 + Y**2 + Z**2) <= 100**2

    metrics_2D_100mm = [None] * n
    metrics_2D_200mm = [None] * n
    for i in range(n):
        m100 = sphere_50  * dis_map_mag_2D[i]; m100 = m100.astype(float); m100[m100 == 0] = np.nan
        m200 = sphere_100 * dis_map_mag_2D[i]; m200 = m200.astype(float); m200[m200 == 0] = np.nan
        metrics_2D_100mm[i] = calc_distortion_metrics(m100)
        metrics_2D_200mm[i] = calc_distortion_metrics(m200)

    m3d_100 = sphere_50  * combined; m3d_100 = m3d_100.astype(float); m3d_100[m3d_100 == 0] = np.nan
    m3d_200 = sphere_100 * combined; m3d_200 = m3d_200.astype(float); m3d_200[m3d_200 == 0] = np.nan

    return {
        'metrics_2D_100mm':    metrics_2D_100mm,
        'metrics_2D_200mm':    metrics_2D_200mm,
        'metrics_3D_est_100mm': calc_distortion_metrics(m3d_100),
        'metrics_3D_est_200mm': calc_distortion_metrics(m3d_200),
        'folders':             folders,
        'qa_image_paths':      qa_image_paths,
    }


def save_csv(results: dict, output_path: str):
    """Save metrics to CSV."""
    folders = results['folders']
    rows = []
    for i, folder in enumerate(folders):
        rows.append({'scan': folder, 'volume': '2D_100mm', **results['metrics_2D_100mm'][i]})
        rows.append({'scan': folder, 'volume': '2D_200mm', **results['metrics_2D_200mm'][i]})
    rows.append({'scan': 'combined', 'volume': '3D_est_100mm', **results['metrics_3D_est_100mm']})
    rows.append({'scan': 'combined', 'volume': '3D_est_200mm', **results['metrics_3D_est_200mm']})

    fieldnames = ['scan', 'volume', 'mean_distortion', 'std_distortion',
                  'max_distortion', 'P95_distortion']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f'Results saved to {output_path}')


def main():
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f'Working directory: {tmpdir}')

        # Authenticate
        xnat_session = get_session(args.url, args.user, args.password)

        scan_labels = [args.scan_axial, args.scan_sagittal, args.scan_coronal]
        scan_dirs = []

        # Download each scan
        for label in scan_labels:
            print(f'Downloading scan: {label}')
            scan_id = get_scan_id(xnat_session, args.url, args.project,
                                  args.subject, args.session, label)
            scan_dir = os.path.join(tmpdir, label)
            os.makedirs(scan_dir, exist_ok=True)
            download_scan(xnat_session, args.url, args.project, args.subject,
                          args.session, scan_id, scan_dir)
            scan_dirs.append(label)

        # Run analysis
        results = run_analysis(tmpdir, scan_dirs, tmpdir)

        # Save and upload CSV
        csv_path = os.path.join(tmpdir, f'distortion_metrics_{args.session}.csv')
        save_csv(results, csv_path)

        print('Uploading results to XNAT...')
        upload_file(xnat_session, args.url, args.project, args.subject,
                    args.session, args.resource_name, csv_path, 'text/csv')

        # Upload QA images
        for img_path in results['qa_image_paths']:
            upload_file(xnat_session, args.url, args.project, args.subject,
                        args.session, args.resource_name, img_path, 'image/png')

    print('Done.')


if __name__ == '__main__':
    main()