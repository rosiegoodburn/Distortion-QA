function [new_grid_mag] = generateDistortionGrid(direc_mr,folder_mr,interp_phantom_grid_gaps)
% Rosie Goodburn, Royal Marsden Hospital

%% Prepare data
% Get list of all DICOM files in folder
dicomFolder = fullfile(direc_mr,folder_mr,'/DICOM');
dicomFiles = dir(fullfile(dicomFolder, '*.dcm'));
fileNames = {dicomFiles.name};
sortedFileNames = natsortfiles(fileNames);

% Preallocate 3D volume matrix
info = dicominfo(fullfile(dicomFolder, sortedFileNames{1}));
vol = zeros(info.Rows, info.Columns, numel(sortedFileNames), 'uint16');

for i = 1:numel(sortedFileNames)
    dicomFilePath_mr_1 = fullfile(dicomFolder, sortedFileNames{i});
    vol(:, :, i) = dicomread(dicomFilePath_mr_1);
end

%% Find grid_pos_mm and orien_mm in scanner frame-of-reference
[orien_mm,grid_pos_mm] = findPhantomPos(vol,info);

%% Crop to phantom grid
mri_cropped = croppingTool(double(vol));

%% Remove bias fields and sum over direction perpendicular to phantom grid
mri_cropped_br = zeros(size(mri_cropped));
for i = 1:size(mri_cropped,1)
    mri_cropped_br(i,:,:) = imflatfield(squeeze(mri_cropped(i,:,:)),15);
end
mri_grid = squeeze(sum(mri_cropped_br,1));

%% Find all points at phantom-grid intersections
pts = gridIntersections_xcorr(mri_grid,8,1,0);

[groupsR,rowIdx,rowY] = groupGridRows(pts,1);  % ±1 px rule
pts_x = pts(:,1).';
pts_y = pts(:,2).';
pts_x(rowIdx == 0) = 0;
pts_y(rowIdx == 0) = 0;
[groupsC,colIdx,colX] = groupGridCols(pts,1);   % ±1 px rule
pts_x(colIdx == 0) = 0;
pts_y(colIdx == 0) = 0;

pts_ = [pts_x; pts_y].';
pts_(any(pts_==0,2),:) = [];

% figure;
% imshow(mri_grid, [], 'InitialMagnification', 'fit');
% hold on
% plot(pts_x, pts_y, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
% hold off

pix_sz = info.PixelSpacing(1); slc_sz = info.SliceThickness;

if size(mri_grid,1) > size(mri_grid,2)
    pts_(:,1) = pts_(:,1) * slc_sz;
    pts_(:,2) = pts_(:,2) * pix_sz;
else
    pts_(:,1) = pts_(:,1) * pix_sz;
    pts_(:,2) = pts_(:,2) * slc_sz;
end

%% Construct 21x21 grid for point locations
step = 14.25;

pts_grid = pts_/step;

% Subtract value of smallest
%temp = mean(pts_grid,2); [~,min_pts_grid_i] = min(temp);
%start = pts_grid(min_pts_grid_i,:);
start = [min(pts_grid(:,1)),min(pts_grid(:,2))];%swapped from pts_grid(1,:)
pts_grid = round(pts_grid + (1-start));

%% Make binary grid
% This will allow orientation
grid_m = zeros(21,21);
for i = 1:size(pts_grid,1)
    %if pts_grid(i,:) > 0
        coords = pts_grid(i,:);
        %vals1 = pts_(i,1); vals2 = pts_(i,2);%???
        grid_m(coords(2),coords(1)) = 1; grid_m(grid_m == 0) = NaN;
    %end
end

%% Build ideal position map for all cross points
[ideal_sideways,ideal_updown] = buildIdealPosMap(orien_mm,grid_pos_mm);

% Detect whether row missing and, if so, remove from grid
no_row_1 = all(isnan(grid_m(1,:))); no_row_2 = all(isnan(grid_m(end,:)));
no_row_3 = all(isnan(grid_m(:,1))); no_row_4  = all(isnan(grid_m(:,end)));
row_missing = no_row_1 || no_row_2 || no_row_3 || no_row_4;

if row_missing
test = ideal_sideways; test(~isnan(test)) = 1; test(isnan(test)) = 0; 
    if sum(sum(test(:,1:9))) == 179
        side_missing = 1;
        ideal_sideways(:,1) = NaN; ideal_updown(1:end,1) = NaN;
    elseif sum(sum(test(1:9,:))) == 179
        side_missing = 2;
        ideal_sideways(1,:) = NaN; ideal_updown(1,:) = NaN;
    elseif sum(sum(test(:,13:end))) == 179
        side_missing = 3;
        ideal_sideways(:,end) = NaN; ideal_updown(:,end) = NaN;
    elseif sum(sum(test(13:end,:))) == 179
        side_missing = 4;
        ideal_sideways(end,:) = NaN; ideal_updown(end,:) = NaN;
    end
end

%% Transform phantom image so that it matches up with orien1_mm
[mri_grid_new] = transformToMatch(mri_grid,grid_m,ideal_sideways);

%% Find all points at phantom-grid intersections
pts = gridIntersections_xcorr(mri_grid_new,8,1,0);

[groupsR,rowIdx,rowY] = groupGridRows(pts,1);  % ±1 px rule
pts_x = pts(:,1).';
pts_y = pts(:,2).';
pts_x(rowIdx == 0) = 0;
pts_y(rowIdx == 0) = 0;
[groupsC,colIdx,colX] = groupGridCols(pts,1);   % ±1 px rule
pts_x(colIdx == 0) = 0;
pts_y(colIdx == 0) = 0;

pts_ = [pts_x; pts_y].';
pts_(any(pts_==0,2),:) = [];

figure;
imshow(mri_grid_new, [], 'InitialMagnification', 1000);
hold on
plot(pts_x, pts_y, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
hold off
pause(0.1);

if size(mri_grid_new,1) > size(mri_grid_new,2)
    pts_(:,1) = pts_(:,1) * slc_sz;
    pts_(:,2) = pts_(:,2) * pix_sz;
else
    pts_(:,1) = pts_(:,1) * pix_sz;
    pts_(:,2) = pts_(:,2) * slc_sz;
end

%% Construct 21x21 grid for point locations
pts_grid = pts_/step;

% Subtract value of smallest
start = [min(pts_grid(1,:)),min(pts_grid(2,:))];
pts_grid = round(pts_grid + (1-start));

%% Make binary grid
% This will allow orientation
grid1 = zeros(21,21); grid2 = zeros(21,21);
grid_m = zeros(21,21);
for i = 1:size(pts_grid,1)
    %if pts_grid(i,:) > 0
        coords = pts_grid(i,:);
        vals1 = pts_(i,1); vals2 = pts_(i,2);%???
        grid_m(coords(2),coords(1)) = 1; grid_m(grid_m == 0) = NaN;
        grid1(coords(2),coords(1)) = vals1;
        grid2(coords(2),coords(1)) = vals2;
    %end
end
grid1 = grid1 .* grid_m; grid2 = grid2 .* grid_m;

%% If phantom-grid row missing, shift grids to match "ideal grid"
if row_missing
test_ideal = ideal_sideways; test_ideal(~isnan(test_ideal)) = 1;
test_ideal(isnan(test_ideal)) = 0;

test_grid = grid1; test_grid(~isnan(test_grid)) = 1;
test_grid(isnan(test_grid)) = 0; 

shifts = [ 0  0;   % no shift
          -1  0;   % up    (row -1)
           1  0;   % down  (row +1)
           0 -1;   % left  (col -1)
           0  1];  % right (col +1)

bestScore = inf;
bestShift = [0 0];

for k = 1:size(shifts,1)
    s = shifts(k,:);
    grids = shift1px_zero_fill(test_grid, s(1), s(2)); % (dRow, dCol)
    score = nnz(xor(grids, test_ideal)); % number of disagreeing pixels
    if score < bestScore
        bestScore = score;
        bestShift = s; % bestShift = [dRow dCol] from binary matching step
    end
end

dRow = bestShift(1); dCol = bestShift(2);

grid1 = shift2d_fill(grid1, dRow, dCol, NaN); % NaN fill recommended
grid2 = shift2d_fill(grid2, dRow, dCol, NaN);

end

%% Interpoalte over phantom-grid gaps
if interp_phantom_grid_gaps == 1
    for i = 1:21
        grid1(:,i) = fillmissing(grid1(:,i),'nearest');
        ideal_sideways(:,i) = fillmissing(ideal_sideways(:,i),'nearest');
        grid2(i,:) = fillmissing(grid2(i,:),'nearest');
        ideal_updown(i,:) = fillmissing(ideal_updown(i,:),'nearest');
    end
    
    if row_missing
        if side_missing == 1
            grid1(:,1) = NaN; grid2(:,1) = NaN;
            ideal_sideways(:,1) = NaN; ideal_updown(:,1) = NaN;
        elseif side_missing == 2
            grid1(1,:) = NaN; grid2(1,:) = NaN;
            ideal_sideways(1,:) = NaN; ideal_updown(1,:) = NaN;
        elseif side_missing == 3
            grid1(:,end) = NaN; grid2(:,end) = NaN;
            ideal_sideways(:,end) = NaN; ideal_updown(:,end) = NaN;
        elseif side_missing == 4
            grid1(end,:) = NaN; grid2(end,:) = NaN;
            ideal_sideways(end,:) = NaN; ideal_updown(end,:) = NaN;
        end
    end
end

%% Correct the positioning of the grids
[grid1,grid2] = correctPosGrids(orien_mm,grid1,grid2,ideal_sideways,ideal_updown);

%% Generate distion maps
distortion_map1 = ideal_sideways - grid1;
distortion_map2 = ideal_updown - grid2;

%% Resample onto 501x501x501 mm volume
new_grid = addDistortionMapToNewGrid(grid_pos_mm,orien_mm,...
    distortion_map1,distortion_map2);

new_grid(isnan(new_grid)) = 0;
new_grid_mag = zeros(size(new_grid,2,3,4));
new_grid_mag(:,:,:) = rssq(new_grid,1);