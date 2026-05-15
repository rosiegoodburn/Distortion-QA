function [orien_mm,grid_pos_mm] = findPhantomPos(volume_mr,info_mr)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

phase_direc = info_mr.InPlanePhaseEncodingDirection;

if strcmp(phase_direc,'ROW')
    phase_direc = 2;
    freq_direc = 1;
else
    phase_direc = 1;
    freq_direc = 2;
end

se2 = strel("sphere",2); se4 = strel("sphere",4);
im_binary = imclose(imopen(imbinarize(double(volume_mr),20),se2),se4);

if sum(im_binary(:)) > numel(im_binary)/2
    im_binary = imclose(imopen(imbinarize(double(volume_mr),220),se2),se4);
end

im_binary = bwareaopen(im_binary,50000);

if freq_direc == 1
    from_top = zeros(size(im_binary,1),size(im_binary,3));
    for i = 1:size(im_binary,1)
        from_top(i,:) = squeeze(sum(im_binary(:,i,:)));
    end
    
    %find the centre of mass (pixel locs) of the knob and main body from looking
    %top-down at phantom
    knob_from_top = from_top;

    knob_from_top(knob_from_top < 0.85 * max(knob_from_top(:))) = 0;
    knob_from_top = bwareaopen(imclose(imbinarize(knob_from_top), se2), 300);
    
    knob_from_top = bwareafilt(knob_from_top,1);

    [r_knob_from_top, c_knob_from_top] = find(knob_from_top == 1);
    rowcolCOM1_knob = [mean(r_knob_from_top), mean(c_knob_from_top)];
    
    main_from_top = imbinarize(from_top);

    [r_main_from_top, c_main_from_top] = find(main_from_top == 1);
    rowcolCOM1_main = [mean(r_main_from_top), mean(c_main_from_top)];
    
    %find the centre of mass (pixel locs) of the knob and main body from looking
    %sideways at the phantom
    crop = round(rowcolCOM1_knob(1));

    %from_side = zeros(crop+1,size(im_binary,3));
    for i = 1:size(im_binary,1)
        from_side(i,:) = squeeze(sum(im_binary(i, crop - 40:crop + 40,:)));
    end
    
    knob_from_side = from_side;

    knob_from_side(knob_from_side > 0.6 * max(knob_from_side(:))) = 0;
    knob_from_side = bwareaopen(imopen(imbinarize(knob_from_side), se4), 300);
    
    [r_knob_from_side, c_knob_from_side] = find(knob_from_side == 1);
    rowcolCOM2_knob = [mean(r_knob_from_side), mean(c_knob_from_side)];

    main_from_side = from_side;

    main_from_side(main_from_side < 0.85 * max(main_from_side(:))) = 0;
    main_from_side = imfill(imbinarize(main_from_side),'holes');
    [r_main_from_side, c_main_from_side] = find(main_from_side == 1);
    rowcolCOM2_main = [mean(r_main_from_side), mean(c_main_from_side)];

    com_knob_pix = [rowcolCOM2_knob(1), rowcolCOM1_knob(1), ...
    mean([rowcolCOM1_knob(2), rowcolCOM2_knob(2)])];

    com_main_pix = [rowcolCOM2_main(1), rowcolCOM1_main(1), ...
    mean([rowcolCOM1_main(2), rowcolCOM2_main(2)])];

elseif freq_direc == 2
    
    from_top = zeros(size(im_binary,1),size(im_binary,3));
    for i = 1:size(im_binary,1)
        from_top(i,:) = squeeze(sum(im_binary(i,:,:)));
    end
    
    %find the centre of mass (pixel locs) of the knob and main body from looking
    %top-down at phantom
    knob_from_top = from_top;

    knob_from_top(knob_from_top < 0.85 * max(knob_from_top(:))) = 0;
    knob_from_top = bwareaopen(imclose(imbinarize(knob_from_top), se2), 300);
    
    [r_knob_from_top, c_knob_from_top] = find(knob_from_top == 1);
    rowcolCOM1_knob = [mean(r_knob_from_top), mean(c_knob_from_top)];
    
    main_from_top = imbinarize(from_top);

    [r_main_from_top, c_main_from_top] = find(main_from_top == 1);
    rowcolCOM1_main = [mean(r_main_from_top), mean(c_main_from_top)];

    %find the centre of mass (pixel locs) of the knob and main body from looking
    %sideways at the phantom
    crop = round(rowcolCOM1_knob(1));

    %from_side = zeros(crop+1,size(im_binary,3));
    for i = 1:size(im_binary,1)
        from_side(i,:) = squeeze(sum(im_binary(crop - 40:crop + 40,i,:)));
    end

    knob_from_side = from_side;
    
    knob_from_side(knob_from_side > 0.6 * max(knob_from_side(:))) = 0;
    knob_from_side = bwareaopen(imopen(imbinarize(knob_from_side), se4), 300);
    
    [r_knob_from_side, c_knob_from_side] = find(knob_from_side == 1);
    rowcolCOM2_knob = [mean(r_knob_from_side), mean(c_knob_from_side)];
    
    main_from_side = from_side;
    main_from_side(main_from_side < 0.85 * max(main_from_side(:))) = 0;
    main_from_side = imfill(imbinarize(main_from_side),'holes');
    [r_main_from_side, c_main_from_side] = find(main_from_side == 1);
    rowcolCOM2_main = [mean(r_main_from_side), mean(c_main_from_side)];
    
    com_knob_pix = [rowcolCOM1_knob(1), rowcolCOM2_knob(1), ...
    mean([rowcolCOM1_knob(2), rowcolCOM2_knob(2)])];
    
    com_main_pix = [rowcolCOM1_main(1), rowcolCOM2_main(1), ...
    mean([rowcolCOM1_main(2), rowcolCOM2_main(2)])];
else
    msg = 'Error: unexpected frequency encoding direction';
    error(msg)
end

orien_dcm = round(info_mr.ImageOrientationPatient);
im_par_pos = info_mr.ImagePositionPatient;
pix_sz = info_mr.PixelSpacing(1);
slc_sz = info_mr.SliceThickness;

vox_sz = [pix_sz,pix_sz,slc_sz];
O2knob_pix = com_knob_pix - vox_sz/2; %O (origin) is [1,1,1]
O2knob_mm = O2knob_pix .* vox_sz;

O2main_pix = com_main_pix - vox_sz/2; %O (origin) is [1,1,1]
O2main_mm = O2main_pix .* vox_sz;

Opos_mm = zeros(1,3);

if all(orien_dcm == [1;0;0;0;1;0])
    
    %[1,1,1] = [-180,-180,171.75](mm), [1.022,1.022,1.5](mm)
    Opos_mm(1) = im_par_pos(1);%?? 1/2 which is which??
    Opos_mm(2) = im_par_pos(2);
    Opos_mm(3) = im_par_pos(3);
    
    knob_pos_mm = zeros(1,3);
    main_pos_mm = zeros(1,3);
    for i = 1:3
        if Opos_mm(i) >= 0
            knob_pos_mm(i) = Opos_mm(i) - O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) - O2main_mm(i);
        else
            knob_pos_mm(i) = Opos_mm(i) + O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) + O2main_mm(i);
        end
    end
    %match to horos %not sure if matching to horos is best - maybe should
    %match to default axial dicom... change this in future!
    knob_pos_mm_temp(1) = knob_pos_mm(2);
    knob_pos_mm_temp(2) = knob_pos_mm(1);
    knob_pos_mm_temp(3) = knob_pos_mm(3);
    knob_pos_mm = knob_pos_mm_temp; clear knob_pos_mm_temp;

    main_pos_mm_temp(1) = main_pos_mm(2);
    main_pos_mm_temp(2) = main_pos_mm(1);
    main_pos_mm_temp(3) = main_pos_mm(3);
    main_pos_mm = main_pos_mm_temp; clear main_pos_mm_temp;

elseif all(orien_dcm == [0;1;0;0;0;-1])

    %[1,1,1] = [-245,180,-170.3](mm), [1.5,1.022,1.022](mm)
    Opos_mm(1) = im_par_pos(3);
    Opos_mm(2) = im_par_pos(2);
    Opos_mm(3) = im_par_pos(1);
    
    knob_pos_mm = zeros(1,3);
    main_pos_mm = zeros(1,3);
    for i = 1:3
        if Opos_mm(i) >= 0
            knob_pos_mm(i) = Opos_mm(i) - O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) - O2main_mm(i);
        else
            knob_pos_mm(i) = Opos_mm(i) + O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) + O2main_mm(i);
        end
    end
    %match to horos
    knob_pos_mm_temp(1) = knob_pos_mm(3);
    knob_pos_mm_temp(2) = knob_pos_mm(2);
    knob_pos_mm_temp(3) = knob_pos_mm(1);
    knob_pos_mm = knob_pos_mm_temp; clear knob_pos_mm_temp;

    main_pos_mm_temp(1) = main_pos_mm(3);
    main_pos_mm_temp(2) = main_pos_mm(2);
    main_pos_mm_temp(3) = main_pos_mm(1);
    main_pos_mm = main_pos_mm_temp; clear main_pos_mm_temp;

elseif all(orien_dcm == [1;0;0;0;0;-1])

    %[1,1,1] = [-185,180,-228](mm), [1.022,1.5,1.022](mm)
    Opos_mm(1) = im_par_pos(3);
    Opos_mm(2) = im_par_pos(1);
    Opos_mm(3) = im_par_pos(2);

    knob_pos_mm = zeros(1,3);
    main_pos_mm = zeros(1,3);
    for i = 1:3
        if Opos_mm(i) >= 0
            knob_pos_mm(i) = Opos_mm(i) - O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) - O2main_mm(i);
        else
            knob_pos_mm(i) = Opos_mm(i) + O2knob_mm(i);
            main_pos_mm(i) = Opos_mm(i) + O2main_mm(i);
        end
    end
    %match to horos
    knob_pos_mm_temp(1) = knob_pos_mm(2);
    knob_pos_mm_temp(2) = knob_pos_mm(3);
    knob_pos_mm_temp(3) = knob_pos_mm(1);
    knob_pos_mm = knob_pos_mm_temp; clear knob_pos_mm_temp;

    main_pos_mm_temp(1) = main_pos_mm(2);
    main_pos_mm_temp(2) = main_pos_mm(3);
    main_pos_mm_temp(3) = main_pos_mm(1);
    main_pos_mm = main_pos_mm_temp; clear main_pos_mm_temp;

else
    msg = 'Error: unexpected orientaion';
    error(msg)
end

%find vector from main body com (mm) to knob com (mm)

main2knob_mm = knob_pos_mm - main_pos_mm;
[~,min_i] = min(abs(main2knob_mm));

%there are 24 possible orientations
orien = zeros(1,3);
for i = 1:3
    orien(i) = sign(main2knob_mm(i));
end

orien_mm = 2*orien;
orien_mm(min_i) = sign(main2knob_mm(min_i));

%calculate com of grid from kob pos + orien (from horos measurement)
if min_i == 2 %horizontal (horos-meas com grid = [1.216, 85.809, -3.750])
    grid2knob_mm = [96.5037, 100.0875, 96.0920];
    grid2knob_mm = grid2knob_mm .* orien;
    grid_pos_mm = knob_pos_mm - grid2knob_mm;
elseif min_i == 3 %vertical1 (horos-meas com grid = [2.360, -57.070, -6.339])
    grid2knob_mm = [96.3468,  98.9717,  99.5187];
    grid2knob_mm = grid2knob_mm .* orien;
    grid_pos_mm = knob_pos_mm - grid2knob_mm;
elseif min_i == 1 %vertical2 (horos-meas com grid = [27.724, -57.075, 6.724])
    grid2knob_mm = [100.2873  99.2393   95.4823];
    grid2knob_mm = grid2knob_mm .* orien;
    grid_pos_mm = knob_pos_mm - grid2knob_mm;
end

%


% %horizontal
% 
% [2,1,2]
% 
% [-2,1,2]
% [2,-1,2]%?? should be impossible
% [2,1,-2]
% 
% [-2,-1,2]%?? should be impossible
% [2,-1,-2]%?? should be impossible
% [-2,1,-2]
% 
% [-2,-1,-2]%?? should be impossible
% 
% %vertical1
% [2,2,1]%this is what the example is
% 
% [-2,2,1]
% [2,-2,1]
% [2,2,-1]
% 
% [-2,-2,1]
% [2,-2,-1]
% [-2,2,-1]
% 
% [-2,-2,-1]
% 
% %vertical2
% [1,2,2]%this is what the example is
% 
% [-1,2,2]
% [1,-2,2]
% [1,2,-2]
% 
% [-1,-2,2]
% [1,-2,-2]
% [-1,2,-2]
% 
% [-1,-2,-2]

end