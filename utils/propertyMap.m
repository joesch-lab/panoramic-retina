function output = propertyMap(idx_map, param)
% calculates mean of parameters for rfs within each bin and returns a
% heatmap across bins

ybins = size(idx_map,1);
xbins = size(idx_map,2);
output = nan(size(idx_map));

for row = 1:ybins
    for col = 1:xbins
        idx = idx_map{row, col}; %DUP
        if numel(idx) < 1; continue; end %no data for this cell
        
        vals = param(idx);

        output(row, col)  = mean(vals, 'all', 'omitnan');
%         output(row, col)  = median(vals, 'all', 'omitnan');
    end
end
end