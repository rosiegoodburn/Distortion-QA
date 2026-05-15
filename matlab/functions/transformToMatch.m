function [mri_grid_new] = transformToMatch(mri_grid,grid,ideal_sideways)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

grid1_1 = grid; %grid1_1(~isnan(grid1_1))=1;
grid1_2 = rot90(grid1_1);
grid1_3 = rot90(grid1_2);
grid1_4 = rot90(grid1_3);

grid1_5 = flip(grid1_1,2);
grid1_6 = rot90(grid1_5);
grid1_7 = rot90(grid1_6);
grid1_8 = rot90(grid1_7);

test = ideal_sideways; test(~isnan(test)) = 1; 

test2 = grid1_1 + test; test3(1) = nansum(test2(:));
test2 = grid1_2 + test; test3(2) = nansum(test2(:));
test2 = grid1_3 + test; test3(3) = nansum(test2(:));
test2 = grid1_4 + test; test3(4) = nansum(test2(:));

test2 = grid1_5 + test; test3(5) = nansum(test2(:));
test2 = grid1_6 + test; test3(6) = nansum(test2(:));
test2 = grid1_7 + test; test3(7) = nansum(test2(:));
test2 = grid1_8 + test; test3(8) = nansum(test2(:));

[~,how_to_transform] = max(test3(:));

if how_to_transform == 1
    mri_grid_new = mri_grid;
elseif how_to_transform == 2
    mri_grid_new = rot90(mri_grid);
elseif how_to_transform == 3
    mri_grid_new = rot90(rot90(mri_grid));
elseif how_to_transform == 4
    mri_grid_new = rot90(rot90(rot90(mri_grid)));
elseif how_to_transform == 5
    mri_grid_new = flip(mri_grid,2);
elseif how_to_transform == 6
    mri_grid_new = rot90(flip(mri_grid,2));
elseif how_to_transform == 7
    mri_grid_new = rot90(rot90(flip(mri_grid,2)));
elseif how_to_transform == 8
    mri_grid_new = rot90(rot90(rot90(flip(mri_grid,2))));
end

end