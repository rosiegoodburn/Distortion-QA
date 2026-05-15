function [ideal_sideways,ideal_updown] = buildIdealPosMap(orien_mm,grid_pos_mm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%construct ideal grid points (!BUT shift to match mean central pnts later)
n = 21;   c = 11;   step = 14.25;
%n = 21;   c = 0;   step = 14.25;
[i, j] = ndgrid(1:n, 1:n);
updown = step * (i - c);   % grid for first dimension
sideways = step * (j - c);   % grid for second dimension

[~,exclude] = min(abs(orien_mm));

dims = [1,2,3]; dims = dims(dims ~= exclude);

if exclude == 2 %horizontal orientation %sideways +ve, updown +ve
elseif exclude == 3 %vertical 1 orientation %sideways +ve, updown -ve
    updown = -1 * updown;
else %%vertical 2 orientation %sideways +ve, updown +ve

end

ideal_grid_sideways_dim_centre = grid_pos_mm(dims(1));
ideal_grid_updown_dim_centre = grid_pos_mm(dims(2));
ideal_sideways = sideways + ideal_grid_sideways_dim_centre;
ideal_updown = updown + ideal_grid_updown_dim_centre;

    %build mask
mask = ones(21,21);
row = [8;9;8;9;10;1;2;3;9;10;11;12;19;20;21;1;2;3;10;11;12;19;20;21;13;14;13;14];
col = [8;8;9;9;9;10;10;10;10;10;10;10;10;10;10;11;11;11;11;11;11;11;11;11;13;13;14;14];
for i = 1:28
    mask(row(i),col(i)) = 0;
end

mask(mask == 0) = NaN;

%there are 8 possible phantom positions for each hor, ver1, ver2
%hor
G{1} = {-2,-1,2}; M{1} = flip(mask,1); %(orien = -2 -1 2)
G{2} = {2,-1,2}; M{2} = rot90(M{1});
G{3} = {2,-1,-2}; M{3} = rot90(M{2});
G{4} = {-2,-1,-2}; M{4} = rot90(M{3});

G{5} = {2,1,2}; M{5} = flip(M{1},2);
G{6} = {2,1,-2}; M{6} = rot90(M{5});
G{7} = {-2,1,-2}; M{7} = rot90(M{6});
G{8} = {-2,1,2}; M{8} = rot90(M{7});

%ver1
G{9} = {-2,-2,-1}; M{9} = flip(mask,1); %(orien = -2 -2 -1)
G{10} = {2,-2,-1}; M{10} = rot90(M{9});
G{11} = {2,2,-1}; M{11} = rot90(M{10});
G{12} = {-2,2,-1}; M{12} = rot90(M{11});

G{13} = {2,-2,1}; M{13} = flip(M{9},2);
G{14} = {2,2,1}; M{14} = rot90(M{13});
G{15} = {-2,2,1}; M{15} = rot90(M{14});
G{16} = {-2,-2,1}; M{16} = rot90(M{15});

%ver2
G{17} = {-1,-2,2}; M{17} = rot90(mask);
G{18} = {-1,2,2}; M{18} = rot90(M{17});
G{19} = {-1,2,-2}; M{19} = rot90(M{18});
G{20} = {-1,-2,-2}; M{20} = rot90(M{19});

G{21} = {1,2,2}; M{21} = flip(M{17},2);
G{22} = {1,2,-2}; M{22} = rot90(M{21});
G{23} = {1,-2,-2}; M{23} = rot90(M{22});
G{24} = {1,-2,2}; M{24} = rot90(M{23});

for i = 1:24
    temp = num2cell(orien_mm);
    test = isequal(temp, G{i});
    if test == 1
        new_mask = M{i};
    end
end

ideal_sideways = ideal_sideways .* new_mask;
ideal_updown = ideal_updown .* new_mask;