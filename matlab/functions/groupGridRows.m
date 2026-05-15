function [groups,rowIdx,rowY] = groupGridRows(pts,tol,minCount)
% pts: N×2 [x y]
% tol: pixels; same-row if |Δy| ≤ tol  (default 1)
% minCount: discard rows with < minCount points (default 8)

if nargin<2 || isempty(tol), tol = 1; end
if nargin<3 || isempty(minCount), minCount = 8; end

y = double(pts(:,2));
g = clusterdata(y,'linkage','single','criterion','distance','cutoff',tol);

u = unique(g);
groups = cell(numel(u),1);
rowY   = zeros(numel(u),1);
for k = 1:numel(u)
    idx = find(g==u(k));
    [~,ord] = sort(pts(idx,1));    % left→right
    groups{k} = idx(ord);
    rowY(k)   = median(y(idx));
end

% order top→bottom
[~,ordR] = sort(rowY);
groups = groups(ordR);
rowY   = rowY(ordR);

% filter by size
keep = cellfun(@numel,groups) >= minCount;
groups = groups(keep);
rowY   = rowY(keep);

% map indices; 0 = discarded
rowIdx = zeros(size(y));
for k = 1:numel(groups), rowIdx(groups{k}) = k; end
end
