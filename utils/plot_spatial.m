function obj = plot_spatial(rf, msg)
rf = squeeze(rf);
colormap(flipud(othercolor('RdBu11')));
obj = imagesc(rf, [-.3 .3]);
xticks([]);
xticklabels([]);
yticks([]);
yticklabels([]);
axis('equal');
axis("square");
axis('xy'); 
xlim([1 size(rf,1)]);
ylim([1 size(rf,2)]);

if nargin > 1
    text(10,10, msg);
end
end

