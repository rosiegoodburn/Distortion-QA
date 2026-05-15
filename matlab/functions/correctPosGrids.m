function [grid1,grid2] = correctPosGrids(orien_mm,grid1,grid2,ideal_sideways,ideal_updown,rowMissing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,test] = min(abs(orien_mm));
if test == 2  %frequency encond == 2 -> horizontal

elseif test == 3 %ver1
    
%figure out how to change all non-nans
    grid2 = -1 * grid2;

    grid2 = grid2 + -1 * (min(grid2(:)) + max(grid2(:)));
elseif test == 1 %ver2
    
end

%figure out where isocentre is...
min_sideways = min(abs(ideal_sideways(:)));
[~,min_sideways_col] = find(abs(ideal_sideways) == min_sideways);
min_updown = min(abs(ideal_updown(:)));
[min_updown_row,~] = find(abs(ideal_updown) == min_updown);

rows = min_updown_row(1)-3:min_updown_row(1)+3;
cols = min_sideways_col(1)-3:min_sideways_col(1)+3;

ideal_sideways_mean = mean(ideal_sideways(rows,cols),"omitnan");
ideal_sideways_mean = mean(ideal_sideways_mean(:),"omitnan");
ideal_updown_mean = mean(ideal_updown(rows,cols),"omitnan");
ideal_updown_mean = mean(ideal_updown_mean(:),"omitnan");

grid1_mean = mean(grid1(rows,cols),"omitnan");
grid1_mean = mean(grid1_mean(:),"omitnan");
grid2_mean = mean(grid2(rows,cols),"omitnan");
grid2_mean = mean(grid2_mean(:),"omitnan");

diff1 = grid1_mean - ideal_sideways_mean;
diff2 = grid2_mean - ideal_updown_mean;

grid1 = grid1 - diff1;
grid2 = grid2 - diff2;

end