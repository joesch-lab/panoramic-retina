% Caluclates the mean of all pixels within rings of different radii
% (thickness 1 px)
function radial_rf = radial_profile(spatial_rf)
n = size(spatial_rf,1);
xc = ceil(n/2);
yc = ceil(n/2);


[xs, ys] = meshgrid(1:n, 1:n);


rmax = n -xc;
radial_rf = zeros(1,rmax);
for r = 1:rmax
    idx = ((xs-xc).^2 + (ys-yc).^2 > r^2) & ((xs-xc).^2 + (ys-yc).^2 < (r+1)^2);
    vals = spatial_rf(idx);
    radial_rf(r) = mean(vals);
end
