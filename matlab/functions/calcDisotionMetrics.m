function [metrics] = calcDisotionMetrics(dis_map_mag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


metrics.mean_distortion = mean(dis_map_mag(:),"omitnan");
metrics.std_distortion = std(dis_map_mag(:),"omitnan");
metrics.max_distortion = max(dis_map_mag(:));
metrics.P95_distortion = prctile(dis_map_mag(:),95);


end