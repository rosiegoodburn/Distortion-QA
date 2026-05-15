% Generates distortion metrics from measurements of the 2D UAL (Uniformity
% and Linearity) phantom.
% 
% 2026, Dr. Rosie Goodburn, Royal Marsden Hospital, UK
%
% Phantom set-up:
%   - Phantom must be completely full of water, including knob compartment
%   - Ensure phantom is not twisted/oblique wrt XYZ
%   - Appropriate setup should ensure that the isocentre lie within the
%     phantom grid
%
% Phantom can be rotated by 90°, flipped, or translated:
%   - Horixontal, perpendicular-to-Y orientations for coronal 2D metrics
%   - Vertical-1, perpendicular-to-Z orientations for axial 2D metrics
%   - Vertical-2, perpendicular-to-X orientations for sagittal 2D metrics

close all
clearvars

scriptFolder = fileparts(mfilename('fullpath'));
targetFolder = fullfile(scriptFolder, 'functions');
dataFolder = '/Users/rgoodburn/Documents/MATLAB/Distortion';
addpath(genpath(targetFolder));

%% Measurements to include:
% direc = fullfile(scriptFolder, 'scans 1');
% folder{1} = '3_fl3d_vibe_PE_RL_AXIAL_horizontal';
% folder{2} = '7_fl3d_vibe_PE_AP_SAGITTAL_vertical1';
% folder{3} = '10_fl3d_vibe_PE_HF_CORONAL_vertical2';

% direc = fullfile(scriptFolder, 'scans 2');
% folder{1} = '4_fl3d_vibe_PE_RL_AXIAL_horizontal';
% folder{2} = '10_fl3d_vibe_PE_AP_SAGITTAL_vertical1';
% folder{3} = '13_fl3d_vibe_PE_HF_CORONAL_vertical2';

% direc = fullfile(scriptFolder, 'scans 3');
% folder{1} = '11_fl3d_vibe_PE_RL_AXIAL_horizontal';
% folder{2} = '19_fl3d_vibe_PE_AP_SAGITTAL_vertical1';
% folder{3} = '22_fl3d_vibe_PE_HF_CORONAL_vertical2';

% direc = fullfile(scriptFolder, 'scans 4');
% folder{1} = '3_fl3d_vibe_PE_RL_AXIAL_horizontal';
% folder{2} = '8_fl3d_vibe_PE_AP_SAGITTAL_vertical1';
% folder{3} = '13_fl3d_vibe_PE_HF_CORONAL_vertical2';

direc = fullfile(dataFolder, 'scans 5');
folder{1} = '3_fl3d_vibe_PE_RL_AXIAL_horizontal';
folder{2} = '7_fl3d_vibe_PE_AP_SAGITTAL_vertical1';
folder{3} = '10_fl3d_vibe_PE_HF_CORONAL_vertical2';

%% Calculate 2D-distortion magnitude map
n = numel(folder);
dis_map_mag_2D = cell(1, n);
for i = 1:n
    dis_map_mag_2D{i} = generateDistortionGrid(direc, folder{i}, 1);
end

%% Estimate 3D-distortion magnitude map
dis_map_mag_3D_est = cell(1, n);
for i = 1:n
    A = dis_map_mag_2D{i};
    B = (A.^2+(A/(2^0.5)).^2).^0.5;
    % Making assumption that e.g., Z^2 = (X^2+Y^2)/2
    dis_map_mag_3D_est{i} = B;
end

%% Combine measurements
% Where overlap, max value taken
combined_dis_map_mag_3D_est = combineArrays(dis_map_mag_3D_est{:});

%% Derive 2D metrics
[X,Y,Z] = meshgrid(-250:250,-250:250,-250:250); 
radius_50 = 50; radius_100 = 100;

sphereVoxels_100 = (X).^2 + (Y).^2 + (Z).^2 <= radius_50.^2;
sphereVoxels_200 = (X).^2 + (Y).^2 + (Z).^2 <= radius_100.^2;

dis_map_mag_2D_100 = cell(1, n); 
dis_map_mag_2D_200 = cell(1, n);

metrics_2D_100mm = cell(1, n);
metrics_2D_200mm = cell(1, n);

for i = 1:n
    dis_map_mag_2D_100{i} = sphereVoxels_100 .* dis_map_mag_2D{i};
    dis_map_mag_2D_200{i} = sphereVoxels_200 .* dis_map_mag_2D{i};
    
    dis_map_mag_2D_100{i}(dis_map_mag_2D_100{i} == 0) = NaN;
    dis_map_mag_2D_200{i}(dis_map_mag_2D_200{i} == 0) = NaN;
    
    metrics_2D_100mm{i} = calcDisotionMetrics(dis_map_mag_2D_100{i});
    metrics_2D_200mm{i} = calcDisotionMetrics(dis_map_mag_2D_200{i});
end

%% Derive 3D metrics
dis_map_mag_3D_est_100 = sphereVoxels_100 .* combined_dis_map_mag_3D_est;
dis_map_mag_3D_est_200 = sphereVoxels_200 .* combined_dis_map_mag_3D_est;

dis_map_mag_3D_est_100(dis_map_mag_3D_est_100 == 0) = NaN;
dis_map_mag_3D_est_200(dis_map_mag_3D_est_200 == 0) = NaN;

metrics_3D_est_100mm = calcDisotionMetrics(dis_map_mag_3D_est_100);
metrics_3D_est_200mm = calcDisotionMetrics(dis_map_mag_3D_est_200);