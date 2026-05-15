function pts = gridIntersections_xcorr(img, L, w, angles)
% img: RGB or grayscale image
% L: half-length of each arm in pixels (≈ 0.6–0.9 of grid spacing)
% w: arm half-width in pixels (≈ grid line half-thickness)
% angles: vector of degrees to handle rotation, e.g. 0 or 0:15:165
% pts: Nx2 [x y] pixel locations

if nargin<3, error('Need L and w'); end
if nargin<4, angles = 0; end

% --- prep ---
if ischar(img) || isstring(img), I = imread(img); else, I = img; end
if size(I,3)==3, I = rgb2gray(I); end
I = im2double(I);
I = 1 - mat2gray(I);            % make dark grid bright

Rmax = zeros(size(I));
for ang = angles
    % --- build normalized cross kernel at this orientation ---
    K = zeros(2*L+1);
    c = L+1;
    K(c, c-L:c+L) = 1;
    K(c-L:c+L, c) = 1;
    if w>0
        K = imdilate(K, strel('diamond', w));
    end
    K = imrotate(K, ang, 'bilinear', 'crop');
    K = imgaussfilt(K, max(1,w));          % smooth arms
    K = K - mean(K(:)); K = K / norm(K(:));% zero-mean, unit-norm

    % --- correlation response ---
    R = imfilter(I, K, 'replicate', 'conv');

    % keep best orientation
    Rmax = max(Rmax, R);
end

SE3 = strel("disk",3); SE6 = strel("disk",6);
mask = imbinarize(I);
mask2 = imerode(mask,SE3);
mask2 = imdilate(mask2,SE6);
mask2 = ~mask2;
% mask2(1:2,:) = 0; mask2(:,1:2) = 0; mask2(end-1:end,:) = 0; mask2(:,end-1:end) = 0;

% --- pick peaks (intersections) ---
Rmax_m = mat2gray(mask2.*Rmax);
M = imextendedmax(Rmax_m, 0.08); % contrast threshold; tune if needed
M = imclearborder(M); %test
M = bwareaopen(M, 1);
M = imclose(M, strel('disk', max(1,round(w/2))));

CC  = bwconncomp(M);
S   = regionprops(CC, Rmax_m, 'WeightedCentroid');
pts = vertcat(S.WeightedCentroid);



end
