function combined_dis_map = combineArrays(varargin)
    % Check inputs
    n = numel(varargin);
    if n == 0
        error('At least one input required.');
    end
    
    arrays = cellfun(@(x) x(:,:,:), varargin, 'UniformOutput', false);
    
    % Stack and average
    combined_dis_map = cat(4, arrays{:});
    combined_dis_map = max(combined_dis_map,[],4,"omitnan");
    
end
