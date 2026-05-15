function [groups,colIdx,colX] = groupGridCols(pts,tol,minCount)
% pts: N×2 [x y]
% tol: pixels; same-col if |Δx| ≤ tol  (default 1)
% minCount: discard columns with < minCount points (default 8)

if nargin<2 || isempty(tol), tol = 1; end
if nargin<3 || isempty(minCount), minCount = 8; end

x = double(pts(:,1));
g = clusterdata(x,'linkage','single','criterion','distance','cutoff',tol);

u = unique(g);
groups = cell(numel(u),1);
colX   = zeros(numel(u),1);
for k = 1:numel(u)
    idx = find(g==u(k));
    [~,ord] = sort(pts(idx,2));    % top→bottom
    groups{k} = idx(ord);
    colX(k)   = median(x(idx));
end

% order left→right
[~,ordC] = sort(colX);
groups = groups(ordC);
colX   = colX(ordC);

% filter by size
keep = cellfun(@numel,groups) >= minCount;
groups = groups(keep);
colX   = colX(keep);

% map indices; 0 = discarded
colIdx = zeros(size(x));
for k = 1:numel(groups), colIdx(groups{k}) = k; end
end
