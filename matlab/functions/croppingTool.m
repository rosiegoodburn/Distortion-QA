function [mri_old_crop2] = croppingTool(mri_old)

around_grid_crop_pix = 3;
top_btm_grid_crop_pix = 3;

if size(mri_old,2) == size(mri_old,3)
    mri_old = permute(mri_old,[3 2 1]);
end
    
%chose top-down view initially
mri_old_sum1 = squeeze(sum(mri_old,1));
mri_old_sum2 = squeeze(sum(mri_old,2));
test1 = max(mri_old_sum1(:)); test2 = max(mri_old_sum2(:));

if test1 < test2
    %top
    mri_old_sum = mri_old_sum1;
    
    crop_loc_rows = gradient(smooth(sum(mri_old_sum,1)));
    [~,crop_loc_rows_1] = max(crop_loc_rows(1:end-3));
    [~,crop_loc_rows_2] = min(crop_loc_rows(1:end-3));
    crop_loc_rows_1 = crop_loc_rows_1 + around_grid_crop_pix;
    crop_loc_rows_2 = crop_loc_rows_2 - around_grid_crop_pix;
    
    crop_loc_cols = gradient(smooth(sum(mri_old_sum,2)));
    [~,crop_loc_cols_1] = max(crop_loc_cols(:));
    [~,crop_loc_cols_2] = min(crop_loc_cols(:));
    crop_loc_cols_1 = crop_loc_cols_1 + around_grid_crop_pix;
    crop_loc_cols_2 = crop_loc_cols_2 - around_grid_crop_pix;
    
    mri_old_crop = mri_old(:,crop_loc_cols_1:crop_loc_cols_2,crop_loc_rows_1:crop_loc_rows_2);
    
    %side
    mri_old_sum2 = squeeze(sum(mri_old_crop,3));
    mask = mri_old_sum2;
    mask(mri_old_sum2 < 0.2 * max(mri_old_sum2(:))) = 0;
    mask(mri_old_sum2 >= 0.2 * max(mri_old_sum2(:))) = 1;
    
    if sum(mask(:))/numel(mask) > 0.3
        mask(mri_old_sum2 < 0.5 * max(mri_old_sum2(:))) = 0;
        mask(mri_old_sum2 >= 0.5 * max(mri_old_sum2(:))) = 1;
    end

    mask = bwareaopen(mask,2000);

    mri_old_sum2 = imflatfield(mri_old_sum2,15,mask);

    crop_loc_freqs1 = smooth(abs(gradient(sum(mri_old_sum2,1))),50);
    crop_loc_freqs2 = smooth(abs(gradient(sum(mri_old_sum2,2))),50);
    
    test1 = abs(crop_loc_freqs1); test1 = max(test1(:));
    test2 = abs(crop_loc_freqs2); test2 = max(test2(:));
    
     if test1 > test2
        %crop_loc_freqs = crop_loc_freqs1;
        %Nfft = 256;
        for a = 1:size(mri_old_sum2,1)
            x(:,a) = smooth(abs(gradient(mri_old_sum2(:,a),2)),50);
        %    [Pxx(a,:),f(a,:)] = pwelch(x(a,:),gausswin(Nfft),Nfft/2,Nfft,1);
        end
        x_ = sum(x,2);
    else
        %crop_loc_freqs = crop_loc_freqs2;
        for a = 1:size(mri_old_sum2,2)
            x(a,:) = smooth(abs(gradient(mri_old_sum2(a,:),1)),50);
        end
        x_ = sum(x,2);
    end
    
    %[pks,locs] = findpeaks(crop_loc_freqs);
    %[~, I] = sort(pks, 'descend');
    
    %j = 0;
    %for i = 1:length(I)      
    %    if I(i) > 2 && I(i) < length(I) - 1 %do not use if I v. close to start/end
    %        j = j + 1;
    %        crop_loc_freqs_test(j) = locs(I(i));
    %    end
    %end

    %test1 = abs(crop_loc_freqs_test(1)-crop_loc_freqs_test(3));
    % test2 = abs(crop_loc_freqs_test(2)-crop_loc_freqs_test(3));
    % if test1 < test2
    %     crop_loc_freqs_(1) = crop_loc_freqs_test(1);
    %     crop_loc_freqs_(2) = crop_loc_freqs_test(3);
    % else
    %     crop_loc_freqs_(1) = crop_loc_freqs_test(2);
    %     crop_loc_freqs_(2) = crop_loc_freqs_test(3);
    % end
    % 
    % crop_loc_freqs_(crop_loc_freqs_ == min(crop_loc_freqs_(:))) = ...
    %     min(crop_loc_freqs_(:)) + plus_n;
    % crop_loc_freqs_(crop_loc_freqs_ == max(crop_loc_freqs_(:))) = ...
    %     max(crop_loc_freqs_(:)) - plus_n;

    crop_loc_freqs_ = find(x_>max(x_)/2);

    mri_old_crop2 = mri_old_crop(min(crop_loc_freqs_) + top_btm_grid_crop_pix : ...
        max(crop_loc_freqs_) - top_btm_grid_crop_pix,:,:);

else
    %top
    mri_old_sum = mri_old_sum2;
    crop_loc_rows = gradient(smooth(sum(mri_old_sum,1)));
    [~,crop_loc_rows_1] = max(crop_loc_rows(1:end-3));
    [~,crop_loc_rows_2] = min(crop_loc_rows(1:end-3));
    crop_loc_rows_1 = crop_loc_rows_1 + around_grid_crop_pix;
    crop_loc_rows_2 = crop_loc_rows_2 - around_grid_crop_pix;
    
    crop_loc_cols = gradient(smooth(sum(mri_old_sum,2)));
    [~,crop_loc_cols_1] = max(crop_loc_cols(:));
    [~,crop_loc_cols_2] = min(crop_loc_cols(:));
    crop_loc_cols_1 = crop_loc_cols_1 + around_grid_crop_pix;
    crop_loc_cols_2 = crop_loc_cols_2 - around_grid_crop_pix;
    
    mri_old_crop = mri_old(crop_loc_cols_1:crop_loc_cols_2,:,crop_loc_rows_1:crop_loc_rows_2);
    
    %side
    mri_old_sum2 = squeeze(sum(mri_old_crop,3));
    mask = mri_old_sum2;
    mask(mri_old_sum2 < 0.2 * max(mri_old_sum2(:))) = 0;
    mask(mri_old_sum2 >= 0.2 * max(mri_old_sum2(:))) = 1;
    
    if sum(mask(:))/numel(mask) > 0.3
        mask(mri_old_sum2 < 0.5 * max(mri_old_sum2(:))) = 0;
        mask(mri_old_sum2 >= 0.5 * max(mri_old_sum2(:))) = 1;
    end

    mask = bwareaopen(mask,2000);

    mri_old_sum2 = imflatfield(mri_old_sum2,15,mask);

    crop_loc_freqs1 = smooth(abs(gradient(sum(mri_old_sum2,1))),50);
    crop_loc_freqs2 = smooth(abs(gradient(sum(mri_old_sum2,2))),50);
    
    test1 = abs(crop_loc_freqs1); test1 = max(test1(:));
    test2 = abs(crop_loc_freqs2); test2 = max(test2(:));
    
    if test1 > test2
        %crop_loc_freqs = crop_loc_freqs1;
        %Nfft = 256;
        for a = 1:size(mri_old_sum2,2)
            x(:,a) = smooth(abs(gradient(mri_old_sum2(:,a),1)),50);
        %    [Pxx(a,:),f(a,:)] = pwelch(x(a,:),gausswin(Nfft),Nfft/2,Nfft,1);
        end
        x_ = sum(x,1);
    else
        %crop_loc_freqs = crop_loc_freqs2;
        for a = 1:size(mri_old_sum2,1)
            x(a,:) = smooth(abs(gradient(mri_old_sum2(a,:),2)),50);
        end
        x_ = sum(x,1);
    end
    
    %[pks,locs] = findpeaks(crop_loc_freqs);
    %[~, I] = sort(pks, 'descend');
    
    %j = 0;
    %for i = 1:length(I)      
    %    if I(i) > 2 && I(i) < length(I) - 1 %do not use if I v. close to start/end
    %        j = j + 1;
    %        crop_loc_freqs_test(j) = locs(I(i));
    %    end
    %end

    %test1 = abs(crop_loc_freqs_test(1)-crop_loc_freqs_test(3));
    % test2 = abs(crop_loc_freqs_test(2)-crop_loc_freqs_test(3));
    % if test1 < test2
    %     crop_loc_freqs_(1) = crop_loc_freqs_test(1);
    %     crop_loc_freqs_(2) = crop_loc_freqs_test(3);
    % else
    %     crop_loc_freqs_(1) = crop_loc_freqs_test(2);
    %     crop_loc_freqs_(2) = crop_loc_freqs_test(3);
    % end
    % 
    % crop_loc_freqs_(crop_loc_freqs_ == min(crop_loc_freqs_(:))) = ...
    %     min(crop_loc_freqs_(:)) + plus_n;
    % crop_loc_freqs_(crop_loc_freqs_ == max(crop_loc_freqs_(:))) = ...
    %     max(crop_loc_freqs_(:)) - plus_n;

    crop_loc_freqs_ = find(x_>max(x_)/2);

    mri_old_crop2 = mri_old_crop(:,min(crop_loc_freqs_) + top_btm_grid_crop_pix : ...
        max(crop_loc_freqs_) - top_btm_grid_crop_pix,:);
end

dims = size(mri_old_crop2);
[~, minDim] = min(dims);
if minDim ~= 1
    order = 1:3;
    order([1, minDim]) = order([minDim, 1]);
    mri_old_crop2 = permute(mri_old_crop2, order);
end

if size(mri_old_crop2,1) > 24
    cut = size(mri_old_crop2,1) - 24;
    if rem(cut, 2) == 0
        cut1 = cut/2; cut2 = cut/2;
    else
        cut1 = floor(cut/2); cut2 = ceil(cut/2);
    end

    mri_old_crop2 = mri_old_crop2(1+cut1:end-cut2,:,:);
end

end