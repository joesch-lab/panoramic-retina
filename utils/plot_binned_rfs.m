function counts = plot_binned_rfs(idx_map, rfs, p)

ybins = size(idx_map, 1);
xbins = size(idx_map, 2);

if nargin < 3;  p = panel(); end
p.pack(ybins, xbins);
p.margin  = 0.2;
counts = [];
for row = 1:ybins
    for col = 1:xbins
        
        idx = idx_map{ybins-row+1, col};
        if numel(idx) < 1; continue; end %no data for this cell
        mean_filter = mean(rfs(:,:,idx),3);
        p(row, col).select();
        %             subplot(nbins, nbins, (row-1)*nbins + col);
        n = numel(idx);
        counts = [counts n];
        if numel(counts) == 1
            plot_spatial(mean_filter, sprintf("n = %d", n));
        else
            plot_spatial(mean_filter, sprintf("%d", n));
        end
    end
end

end