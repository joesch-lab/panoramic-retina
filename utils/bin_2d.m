
function idx_map = bin_2d(centroids, varargin)

parser = inputParser;
parser.addRequired('centroids');
parser.addParameter('binsize', 100); %size of bins in microns
parser.addParameter('mincount', 5); %minimum neurons per bin
parser.addParameter('lims', [-1500 1500]);
parser.parse(centroids, varargin{:});
mincount = parser.Results.mincount;
binsize = parser.Results.binsize;
lims = parser.Results.lims;

xmin = lims(1); ymin = xmin;
xmax = lims(2); ymax = xmax;
xedg = xmin:binsize:xmax;
yedg = ymin:binsize:ymax;
xbins = numel(xedg) - 1;
ybins = numel(yedg) - 1;

for row = 1:ybins
    for col = 1:xbins
        x1 = xedg(col); x2 = xedg(col+1);
        y1 = yedg(row); y2 = yedg(row+1);
        
        xidx = centroids(:,1) > x1 & centroids(:,1) < x2;
        yidx = centroids(:,2) > y1 & centroids(:,2) < y2;
        
        idx = find(yidx & xidx);
        if nnz(idx) < mincount %too few cells in this bin, skip
            idx_map{row,col} = [];
        else
            idx_map{row, col} = idx;
        end
    end
end
end