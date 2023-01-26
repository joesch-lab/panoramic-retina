function pts = plot_trend_by_group(centroids, data, group_id)
nb = 6;
pts = [];
op.nBins = nb;
op.LineWidth = 2;
op.DisplayName = sprintf("All Neurons, n = %d", numel(data));
plot_conditional_linear(centroids(:,2), data, op);
hold on;
groups = unique(group_id);
colors = othercolor('Dark28', numel(groups));
for j = 1:numel(groups)
    i = groups(j);
    idx = group_id == i;
    opts.plotSEM = false;
    opts.meancolor = colors(numel(groups)+ 1-j,:);
    opts.nBins = nb;
    opts.DisplayName = sprintf("Retina %d, n = %d", i, nnz(idx));
    temp = plot_conditional_linear(centroids(idx, 2), data(idx), opts);
    pts = [pts; temp];
    xlim([-1500 1500]);
    xticks([-1500 0 1500]);
end
legend('Location', 'northeast');
xlabel("<- Dorsal Ventral ->");

    function pts = plot_conditional_linear(x, y, opts)
        %% Conditional probability of y at different bins of x
        % plots mean +- SEM
        
        if nargin < 3  || ~isfield(opts, 'nBins')
            nBins = 10;
        else
            nBins = opts.nBins;
        end
        if nargin < 3 || ~isfield(opts, 'plotSEM')
            plotSEM = true;
        else
            plotSEM = opts.plotSEM;
        end
        if nargin < 3 || ~isfield(opts, 'LineWidth')
            linewidth = 1;
        else
            linewidth = opts.LineWidth;
        end
        if nargin < 3 || ~isfield(opts, 'meancolor')
            meancolor = [0 0 0];
        else
            meancolor = opts.meancolor;
        end
        if nargin < 3 || ~isfield(opts, 'DisplayName')
            dispname = "";
        else
            dispname = opts.DisplayName;
        end
        if sum(isnan(y)) > 0 || nnz(y==Inf) > 0%nanproof the function
            keep = ~isnan(y) & ~(y==Inf);
            y = y(keep);
            x = x(keep,:);
        end

        binsV = linspace(min(x), max(x), nBins + 1);
        avgV = zeros(nBins, 1);
        semV = zeros(nBins, 1);
        for n = 1:nBins
            values = y(x >= binsV(n) & x < binsV(n+1));
            avgV(n) = mean(values);
            semV(n) = std(values) / sqrt(numel(values));
        end
        
        x_vals = (binsV(1:end-1) + diff(binsV)/2)';
        plot(x_vals, avgV, 'Color', meancolor, 'DisplayName', dispname, 'LineWidth', linewidth);
        y_top = avgV + semV;
        y_bot = avgV - semV;
        hold on;
        if plotSEM
            hold on;
            patch([x_vals; flipud(x_vals)], [y_top; flipud(y_bot)], 'b', 'FaceAlpha', 0.2, 'DisplayName', 'SEM', 'LineStyle', 'none');
        end
        pts = [x_vals(:) avgV(:)];
    end
end