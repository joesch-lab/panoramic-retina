
function norm_filter = normalize_filter(filter, filter_opt)
if nargin < 2; filter_opt = []; end 
% negative_times = filter_opt.latencies <= 0;
filter = squeeze(filter);

if ndims(filter) == 2
    dims = [1,2];
    w = size(filter, 1);
    h = size(filter, 2);
    [x,y] = meshgrid(1:w,1:h);
    border = x < .1*w | x > .9*w | y < .1*h | y > .9*h;
    baseline = mean(filter(border), 'all');
    norm_filter = filter - baseline;
elseif ndims(filter) == 3 % x*y*time
    negative_times = 1:5;
    dims = [1,2,3];
    norm_filter = filter - mean(filter(:,:,negative_times), dims);
elseif ndims(filter) == 4 %neuron*x*y*time
    negative_times = 1:5;
    dims = [2,3,4];
    norm_filter = filter - mean(filter(:,:,:,negative_times), dims);
end

norm_filter = norm_filter./max(abs(norm_filter), [], dims);
end