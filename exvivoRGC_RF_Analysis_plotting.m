%Tested on Matlab 2019b

%requires:
% panel (from Matlab file exchange)
% othercolor (from Matlab file exchange)


load("./exvivo_RFs.mat");
addpath("./utils");


%%
fprintf("Flipping ON rfs to OFF\n");
spatial_unsigned_all = flip_on_rfs(spatial_all);
clear spatial_all; %save RAM


%% Fig 3a. Binned Spatial RFs from one retina (retina 5)

f = figure;
idx_ret = find(ret_number == 5);
idx_map = bin_2d(centroids(idx_ret,:), 'binsize', 300, 'mincount',5);
counts = plot_binned_rfs(idx_map, spatial_unsigned_all(:,:,idx_ret));
sgtitle(sprintf("Fig 3a: Cell count = %.1f ± %.1f", mean(counts), std(counts)));
f.Position = [0 0 750 750]; %ensure figure is square for proper alignment

%% Extended Data Fig 5. Binned Spatial RFs from all sopsin retinas

f = figure;
idx_ret = find(sopsin_gt);
idx_map = bin_2d(centroids(idx_ret,:), 'binsize', 300, 'mincount',5);
counts = plot_binned_rfs(idx_map, spatial_unsigned_all(:,:,idx_ret));
sgtitle(sprintf("Extended Fig 5: Cell count = %.1f ± %.1f", mean(counts), std(counts)));
f.Position = [0 0 750 750]; %ensure figure is square for proper alignment

%% Fig 3b. Mean RFs from below, at and above horizon

figure;
cond = sopsin_gt; %only use sopsin retinas
as = [-1500, -900, 500]; bs = [-1000, -700, 600]; %bin limits in um
for i = 1:3
    subplot(2,3,3+i);
    a = as(i); b = bs(i);
    idx = centroids(:,2) > a & centroids(:,2) < b & cond;
    spatial = squeeze(mean(spatial_unsigned_all(:,:,idx), 3));
    plot_spatial(spatial, sprintf("N = %d", nnz(idx)));
    title(sprintf("[%d, %d] µm", a, b));
end
sgtitle("Fig 3b");



%% Fig 3c. Radial RF profiles (same bins as above)

fprintf("computing radial profiles of all rfs, ~2mins\n");
radial_all = zeros(size(spatial_unsigned_all,3), 50);
for i = 1:size(spatial_unsigned_all,3)
    radial_all(i,:) = radial_profile(spatial_unsigned_all(:,:,i));
end
%%
for i = 1:3
    a = as(i); b = bs(i);
    idx = centroids(:,2) > a & centroids(:,2) < b & cond;
    radials = radial_all(idx,:);
    avg_rad(i,:) = mean(radials,1);
    caption{i} = sprintf("[%d, %d] µm, n = %d", a, b, nnz(idx));
end

figure;
subplot(2,1,2);
plot(linspace(0,500,size(radials,2)),avg_rad', 'LineWidth', 3);
yline(0, 'LineStyle', '--');
legend(caption,'Location', 'southeastoutside');
xlabel("Radius (um)");
ylabel("RF Profile Weight");
ylim([-.1 .1]);
yticks([-.1 0 .1]);
xticks([0 500]);
xlim([0 500]);
colorv = [144 25 28]/255;
colorh = [184 72 38]/255;
colord = [231 138 36]/255;
colororder([colord;colorh;colorv]);
sgtitle("Fig 3c");

%% Fig 3d-i RF Properties binned in 1D or 2D
% and extended Data Fig. 7

clear cond;
cond{1} = find(sopsin_gt);
cond{2} = find(~sopsin_gt);

for c = 1:2    
    idx = cond{c};
    idx_map = bin_2d(centroids(idx,:), 'binsize', 100);
    colors = {flipud(othercolor('PRGn7')), flipud(othercolor('PRGn7')), othercolor('BrBG11')};
    titles = ["Relative Surround Strength", "Center Size", "Vertical Surround Asymmetry"];
    params{1} = scRatio(idx);
    params{2} = centSize(idx);
    params{3} = rfAsymm(idx);
    
    figure;
    if c == 1; sgtitle("Fig 3d-i"); else sgtitle("Extended Fig. 7"); end
    for i = 1:3
        % 2D Bins
        subplot(2,3,i);
        output = propertyMap(idx_map, params{i});
        nz = output(output~=0);
        mini = quantile(nz(:), .05); %remove extremes from dynamic range
        maxi = quantile(nz(:), .95);
        title(titles(i));
        imAlpha = ones(size(output));
        imAlpha(isnan(output)) = 0;
        hold on;
        if i~=3
            imagesc(output, 'AlphaData', imAlpha, [mini, maxi]);
        else
            imagesc(output, 'AlphaData', imAlpha, [-1,1]);
        end
        xticks([]); yticks([]);
        colorbar;
        ax = imgca();
        axis image; axis xy;
        colormap(ax, colors{i});
        mid = size(output,1)/2;
        plot(mid,mid, '+');
        
        % 1D Bins (per retina)
        subplot(2,3,3+i);
        pts = plot_trend_by_group(centroids(idx,:), params{i}, ret_number(idx));
        ylabel(titles(i));
        if i~=3; legend('hide'); end
        dors = pts(:,1) < 0;
        [~,pval] = kstest2(pts(dors,2), pts(~dors,2));
        title(sprintf("Pvalue: %e", pval));
    end
end


%% Fig. 3j (Example surround orientations from retina 9)
figure;

colors = [colorv;colorh;colord];

for i = 1:3
    subplot(3,1,i);
    
    orts = [gprops_all(idxs{i}).orientation_centers_of_mass_data_masks];
    alpha = wrapTo180(sopsin_dorsal(9)-90);
    orts = wrapToPi(-deg2rad(alpha) - orts);
    
    histogram(orts, 20,'Normalization', 'probability', 'FaceColor', colors(i,:), 'EdgeAlpha', 0);
    ylim([0 .5]);
    hold on;
    xlim([-pi pi]);
    xticks(linspace(-pi,pi, 5));
    xticklabels(["T", "D", "N", "V", "T"]);
end

sgtitle("Fig. 3j");


%% Fig. 3k (Surround asymmetry in sinusoidal projection of visual space)

alpha0 = deg2rad(22); %final coordinates of optical axis = resting position of mouse eye
theta0 = deg2rad(64);

cents = polar_to_visual(polar_coords, alpha0, theta0);
minlim = min(cents, [], 'all'); maxlim = max(cents, [], 'all');
lims = [minlim maxlim];
binsize = (maxlim - minlim)/40; %divide into 40 bins
idx_map = bin_2d(cents, 'binsize', binsize, 'lims', lims, 'mincount', 0); 
map = propertyMap(idx_map, rfAsymm);

figure;
colormap(othercolor('BrBG11'));
imAlpha = ones(size(map));
imAlpha(isnan(map)) = 0;
image('CData', map, 'XData', lims, 'YData', lims,'CDataMapping','scaled', 'AlphaData', imAlpha);
caxis([-1 1]);
colorbar;
hold on;
plotGrid();
axis equal

% Also plot the rim of a retina
phiRim = zeros(1,100)'; %latitude = 0 defines the equator in retinal coords
lambdaRim = linspace(-pi,pi,100)'; % make a line passing through all longitudes
cent_rim = polar_to_visual([lambdaRim phiRim], alpha0, theta0);
x = cent_rim(:,1); y = cent_rim(:,2);
kink = find(abs(diff(x)) > 2); 
x = circshift(x,-kink);
y = circshift(y, -kink);
plot(x,y, 'LineWidth', 2, 'Color', 'k');

title("Surround Asymmetry in Visual Coordinates");

xticks([-pi 0 pi]);
xticklabels({'-180', 'Front', '180'});
yticks([-pi/2 0 pi/2]);
yticklabels({'ground', 'horizon', 'sky'});
xlabel("Azimuth");
ylabel("Elevation");


%% Fig 4a

k = numel(unique(group_idx));

temps = temporal_all./max(abs(temporal_all), [], 2); %normalize


nCellsPerClust = zeros(k,1);

for i = 1:k
    nCellsPerClust(i) = nnz(group_idx == i);
end
inds = zeros(length(group_idx), 1);
idxPic = zeros(length(group_idx), 1);
ccmap = jet(k);
cnt = 0;
for cluster_id = 1:k
    cnt1 = cnt + sum(group_idx == cluster_id);
    inds(cnt+1:cnt1) = find(group_idx == cluster_id);
    idxPic(cnt+1:cnt1) = cluster_id;     
    cnt = cnt1;
end

figure;
ax_cluster = subplot(1, 12, 1);
imagesc(idxPic);
ylabel('Cluster ID'); 
xlabel('');
set(gca, 'YTickLabels', {''}, 'YTick', []);
colormap(ax_cluster, jet)

ax_data = subplot(1, 12, 2:12);
imagesc(temps(inds,:), [-1.5,1.5]);
colormap(flipud(othercolor('RdBu11')));
xticks([1 size(temps,2)]);
xticklabels({0 num2str(25*size(temps,2))});
title('Temporal RFs'); 
xlabel('Time (ms)');
sgtitle("Fig 4a");

figure;
for id = 1:k
    subplot(5,2,id);
    idx = group_idx == id;
    traces = temps(idx,:);
    avg = mean(traces, 1)';
    stds = std(traces, 0, 1)';
    hold on;
    plot(avg-stds, 'Color', 'k', 'LineWidth', 1); 
    plot(avg+stds, 'Color', 'k', 'LineWidth', 1);
    plot(avg, 'LineWidth', 2, 'Color', "k");
    title(sprintf("cluster %d, n = %d", id,nnz(group_idx == id)));
    ylim([-1.2 1.2]);
    xlim([1, size(temps, 2)]);
    xticks([1 size(temps,2)]);
    xticklabels({0 num2str(25*size(temps,2))});
end
sgtitle("Fig 4a");




%% Fig 4c-e 1D Bins of properties, within each cluster

idx = find(sopsin_gt);
titles = ["Relative Surround Strength", "Center Size", "Vertical Surround Asymmetry"];
params{1} = scRatio(idx);
params{2} = centSize(idx);
params{3} = rfAsymm(idx);

figure;
for i = 1:3    
    subplot(1,3,i);
    pts = plot_trend_by_group(centroids(idx,:), params{i}, group_idx(idx));
    ylabel(titles(i));
    if i~=3; legend('hide'); end
    dors = pts(:,1) < 0;
    [~,pval] = kstest2(pts(dors,2), pts(~dors,2));
    title(sprintf("Pvalue: %e", pval));
end





%% Helper functions

function unsigned = flip_on_rfs(signed)
unsigned = signed;
mid = ceil(size(signed,1)/2);
onrfs = signed(mid,mid,:) > 0;
unsigned(:,:,onrfs) = -signed(:,:,onrfs);
end
