function [new_grid] = addDistortionMapToNewGrid(grid_pos_mm,orien_mm,...
    distortion_map1,distortion_map2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

new_grid = zeros(3,501,501,501);
[Xq,Yq,Zq] = meshgrid(-250:250,-250:250,-250:250);


[~,test] = min(abs(orien_mm));
if test == 2  %frequency encond == 2 (Y) -> horizontal orientation
    
    %MATALB: First dimension = rows (Y-axis). (updown)
            %Second dimension = columns (X-axis). (sideways)
            %Third dimension = slices or planes (Z-axis).
    
    %The distortion map corresponds to:
    % sideways = RL +ve, X scanner dimension -> 2nd dimension in map??
    % updown = IS +ve, Z scanner dimension -> 1st dimention in map??
    % slices = PA -ve (?), Y scanner dimension
    
    distortion_map_ = NaN(3,3,size(distortion_map1,1)+2,size(distortion_map1,2)+2);
    for i = 1:3
        distortion_map_(1,i,2:end-1,2:end-1) = distortion_map1; %check this distortion direciton
        distortion_map_(3,i,2:end-1,2:end-1) = distortion_map2; %check this distortion direciton        
    end

    distortion_map_ = permute(distortion_map_,[1 2 4 3]);

    V1 = squeeze(distortion_map_(1,:,:,:));
    V3 = squeeze(distortion_map_(3,:,:,:));
    %sign_orien = sign(orien_mm);
    
    X_start = grid_pos_mm(1)-14.25*11;
    X_end = grid_pos_mm(1)+14.25*11;
    X_step = 14.25 * sign(X_end - X_start);

    Y_start = grid_pos_mm(2)-14.25;
    Y_end = grid_pos_mm(2)+14.25;
    Y_step = 14.25 * sign(Y_end - Y_start);

    Z_start = grid_pos_mm(3)-14.25*11;
    Z_end = grid_pos_mm(3)+14.25*11;
    Z_step = 14.25 * sign(Z_end - Z_start);

    [X,Y,Z] = meshgrid(X_start:X_step:X_end, Y_start:Y_step:Y_end, Z_start:Z_step:Z_end);
    
    V1q = interp3(X,Y,Z,V1,Xq,Yq,Zq,'nearest');
    V3q = interp3(X,Y,Z,V3,Xq,Yq,Zq,'nearest');
    
    new_grid(1,:,:,:) = V1q;
    new_grid(3,:,:,:) = V3q;

elseif test == 3 %frequency encond == 3 (Z) -> vertical-1 orientation

    %MATALB: First dimension = rows (Y-axis). (updown)
            %Second dimension = columns (X-axis). (sideways)
            %Third dimension = slices or planes (Z-axis). (inout)
    
    %The distortion map corresponds to:
    % sideways = RL +ve, X scanner dimension -> 2nd dimension in map??
    % updown = PA -ve, Y scanner dimension -> 1st dimention in map??
    % slices = IS -ve (?), Z scanner dimension
    
    distortion_map_ = NaN(3,3,size(distortion_map1,1)+2,size(distortion_map1,2)+2);
    for i = 1:3
        distortion_map_(1,i,2:end-1,2:end-1) = distortion_map1; %check this distortion direciton
        distortion_map_(2,i,2:end-1,2:end-1) = distortion_map2; %check this distortion direciton        
    end

    distortion_map_ = permute(distortion_map_,[1 3 4 2]);
    
    V1 = squeeze(distortion_map_(1,:,:,:));
    V2 = squeeze(distortion_map_(2,:,:,:));

    X_start = grid_pos_mm(1)-14.25*11;
    X_end = grid_pos_mm(1)+14.25*11;
    X_step = 14.25 * sign(X_end - X_start);

    Y_start = grid_pos_mm(2)+14.25*11;
    Y_end =  grid_pos_mm(2)-14.25*11;
    Y_step = 14.25 * sign(Y_end - Y_start);

    Z_start = grid_pos_mm(3)+14.25;
    Z_end = grid_pos_mm(3)-14.25;
    Z_step = 14.25 * sign(Z_end - Z_start);

    [X,Y,Z] = meshgrid(X_start:X_step:X_end, Y_start:Y_step:Y_end, Z_start:Z_step:Z_end);
    
    V1q = interp3(X,Y,Z,V1,Xq,Yq,Zq,'nearest');
    V2q = interp3(X,Y,Z,V2,Xq,Yq,Zq,'nearest');
    
    new_grid(1,:,:,:) = V1q;
    new_grid(2,:,:,:) = V2q;

    % %grid_pos_mm(3) = orien_mm(3) * grid_pos_mm(3);
    % grid_pos = uint16(round(grid_pos_mm) + size_new_grid/2 + 1);
    % 
    % grid_pos(1) = grid_pos(1) - 150;
    % grid_pos(2) = grid_pos(2) - 150;
    % grid_pos(3) = grid_pos(3) - 16;
    % 
    % new_grid(1,grid_pos(2):grid_pos(2)+298,grid_pos(1):grid_pos(1)+298,...
    %     grid_pos(3):grid_pos(3)+30) = distortion_map_(1,:,:,:);
    % new_grid(2,grid_pos(2):grid_pos(2)+298,grid_pos(1):grid_pos(1)+298,...
    %     grid_pos(3):grid_pos(3)+30) = distortion_map_(2,:,:,:);
    % 
    % new_grid_for_viewing = flip(imrotate3(squeeze(new_grid(1,:,:,:)),-90,[0 0 1]),2);

elseif test == 1 %frequency encond == 1 (X) -> vertical-2 orientation
    
    %MATALB: First dimension = rows (Y-axis). (updown)
            %Second dimension = columns (X-axis). (sideways)
            %Third dimension = slices or planes (Z-axis). (inout)
    
    %The distortion map corresponds to: %fix!!!
    % sideways = AP +ve, Y scanner dimension -> 2nd dimension in map??
    % updown = IS +ve, Z scanner dimension -> 1st dimention in map??
    % slices = RL -ve (?), X scanner dimension 

    distortion_map_ = NaN(3,3,size(distortion_map1,1)+2,size(distortion_map1,2)+2);
    for i = 1:3
        distortion_map_(2,i,2:end-1,2:end-1) = distortion_map1; %check this distortion direciton
        distortion_map_(3,i,2:end-1,2:end-1) = distortion_map2; %check this distortion direciton        
    end

    distortion_map_ = permute(distortion_map_,[1 4 2 3]);
        
    V2 = squeeze(distortion_map_(2,:,:,:));
    V3 = squeeze(distortion_map_(3,:,:,:));

    X_start = grid_pos_mm(1)-14.25;
    X_end = grid_pos_mm(1)+14.25;
    X_step = 14.25 * sign(X_end - X_start);

    Y_start = grid_pos_mm(2)-14.25*11;
    Y_end =  grid_pos_mm(2)+14.25*11;
    Y_step = 14.25 * sign(Y_end - Y_start);

    Z_start = grid_pos_mm(3)-14.25*11;
    Z_end = grid_pos_mm(3)+14.25*11;
    Z_step = 14.25 * sign(Z_end - Z_start);

    [X,Y,Z] = meshgrid(X_start:X_step:X_end, Y_start:Y_step:Y_end, Z_start:Z_step:Z_end);
    
    V2q = interp3(X,Y,Z,V2,Xq,Yq,Zq,'nearest');
    V3q = interp3(X,Y,Z,V3,Xq,Yq,Zq,'nearest');
    
    new_grid(2,:,:,:) = V2q;
    new_grid(3,:,:,:) = V3q;


    % grid_pos_mm(1) = orien_mm(1) * grid_pos_mm(1);
    % grid_pos = uint16(round(grid_pos_mm) + floor(size_new_grid/2) + 1);
    % 
    % grid_pos(1) = grid_pos(1) - 16;
    % grid_pos(2) = grid_pos(2) - 150;
    % grid_pos(3) = grid_pos(3) - 150;
    % 
    % new_grid(2,grid_pos(2):grid_pos(2)+298,grid_pos(1):grid_pos(1)+30,...
    %     grid_pos(3):grid_pos(3)+298) = distortion_map_(2,:,:,:);
    % new_grid(3,grid_pos(2):grid_pos(2)+298,grid_pos(1):grid_pos(1)+30,...
    %     grid_pos(3):grid_pos(3)+298) = distortion_map_(3,:,:,:);

end

% % INPUTS you must supply
% A    = distortion_map_;                % size 31x299x299   (rows, cols, slices)
% R    = orien_mm;           % 3x3, columns = unit axes of source (x,y,z) in world
% com  = grid_pos_mm(:);                % 3x1 world coords of the center you want at (0,0,0)
% voxS = [1 1 1];                % source voxel sizes in world units: [row, col, slice]
% voxT = 1;                         % desired target voxel size in same units (scalar)
% 
% % Sizes and centers (MATLAB order: [row col slice])
% SzS = size(A);                    % [31 299 299]
% cS  = (SzS+1)/2;                  % source index center (1-based), e.g. [16 150 150]
% SzT = [501 501 501];
% cT  = (SzT+1)/2;                  % [251 251 251]
% 
% % Build linear map  (indices -> world -> target indices)
% % Step 1: source indices -> world (center, scale, rotate)
% Sscale = diag(voxS);              % scale rows, cols, slices into world units
% Mlin   = (1/voxT) * R * Sscale;   % world -> divide by target spacing
% % Step 2: full affine to target indices (add centers and COM shift)
% tvec = cT(:) - Mlin*cS(:) - (1/voxT)*com(:);  % 3x1
% 
% % Assemble 4x4 for affine3d (affine3d expects post-multiply of row vectors)
% T = eye(4);
% T(1:3,1:3) = Mlin.';
% T(4,1:3)   = tvec.';   % translation in last row for row-vector convention
% 
% tform = affine3d(T);
% 
% % Output spatial ref (axis-aligned; origin 0 at voxel [251 251 251])
% % imref3d only needs voxel spacing; world limits are implied by cT
% RB = imref3d(SzT, voxT, voxT, voxT);
% 
% % Resample
% A501 = imwarp(A, tform, 'OutputView', RB, 'InterpolationMethod','linear', 'FillValues', 0);



end