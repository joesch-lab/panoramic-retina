
% Plotting code for Ext. Fig. 4
% assumes you have already run RF_gaussian_mixture.m and gaussian_mixture_RF_properties.m
% and the RFprops are saved in RFprops.mat


gprops = load('RFprops.mat') 

cnt = 1;
z=randsample(numel(gprops), 100);

%Scroll through a random sample 
for i = z'
    clf;
    subplot(2,3,[1,4]);
    gprop = gprops(i);
    im = gprop.RFaroundpeak;
    imagesc(im, [-.5 .5]);
    xticks([]); yticks([]);
    colormap(flipud(othercolor('RdBu11')));
    title("RF");
    axis off;
    axis image;
    
    subplot(2,3,[2 5]);
    model = gprop.model;
    maxamp=max(max(model(:)),abs(min(model(:))));
    modeli=model./maxamp;      % in the range[-1,1]
    imagesc(modeli, [-.5 .5]);
    title("Model");
    xticks([]); yticks([]);
    axis image;
    axis off;
    colormap(flipud(othercolor('RdBu11')));

    %plot ellipses
    hold on;
    center = gprop.center;
    surround = gprop.surround;
    if isempty(surround); continue; end;
    
    if center(1)>0
        
        plot_ellipse(center(2),center(3),center(4)*2,center(5)*2,-center(6),[1,0,0],0.0);
        if surround(1)~=0
            hold on;
            plot_ellipse(surround(2),surround(3),surround(4)*2,surround(5)*2,-surround(6),[0,0,1],0.0);
        end
    else
        plot_ellipse(center(2),center(3),center(4)*2,center(5)*2,-center(6),[0,0,1],0.0);
        if surround(1)~=0
            hold on;
            plot_ellipse(surround(2),surround(3),surround(4)*2,surround(5)*2,-surround(6),[1,0,0],0.0);
        end
    end
    

    subplot(2,3,3);
    plot(gprop.center_av_value_data);
    ylim([-1 1]);
    xticks([40]);
    yticks([-1 0 1]);
    title("Center Dynamics");
            
    
    subplot(2,3,6);
    plot(gprop.surround_av_value_data);
    ylim([-.2 .2]);
    xticks([40]);
    yticks([-.2 0 .2]);
    title("Surround Dynamcis");

    rf = gprop.RFaroundpeak;
    smask = logical(gprop.surround_mask);
    cmask = logical(gprop.center_mask);

    centSize = nnz(cmask);
    surrSize = nnz(smask);
    centStrength = abs(sum(rf(cmask), 'all'));
    surrStrength = abs(sum(rf(smask), 'all'));
    scRatio = surrStrength./centStrength;
    if scRatio > 10
        scRatio = 10; %surround 10x stronger than center, unrealistic
    end

    rf(~smask) = 0;
    rf(cmask) = 0;
    upper = rf(1:101<51,:);
    lower = rf(1:101>51,:);

    u = abs(sum(upper, 'all'));
    l = abs(sum(lower, 'all'));
    rfAsymm = (l-u)/(u+l);

    
    tit = [];
    tit = [tit sprintf("SNR: %f, R2: %f", gprop.peak2noise_log, gprop.R2_cs)];
    tit = [tit sprintf(" centSize: %f, scRatio: %f, rfasymm: %f",  centSize, scRatio, rfAsymm)];
    tit = [tit sprintf("d: %f, ort: %f", gprop.distance_btw_centers_of_mass, gprop.orientation_centers_of_mass)];
    sgtitle(tit);
    
    str = input("Press enter to continue: ", 's');
end
