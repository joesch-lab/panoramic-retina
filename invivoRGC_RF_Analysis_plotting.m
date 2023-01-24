%% %%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elmin=-35; %lowest elevation (deg)
elmax=35; %highest elevation (deg)
nbins_low=16; %number of bins for coarse binning
nbins_high=320; %number of bins for fine binning
mean_filt_width=5; % gaussian filter width for visualization (fine binning)  (deg)
statbins_el=[-inf 0 5 inf]; % bins for ks test comparison elevation (1st vs. 3rd) (deg)
statbins_az=[-inf 60 65 inf];% bins for ks test comparison azimuth (1st vs. 3rd) (deg)
minRGCparam=18;


%% %%%%%%%%%%%%%%%%%%%%%%%%% LOADING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\RGC_in_vivo_RFs_selected.mat');
figdir='D:\RGCfigs\';mkdir(figdir)
%% %%%%%%%%%%%%%%%%%%%%%%%%% COMPUTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elevation and azimuth bins
bins_el_low = linspace(elmin,elmax,nbins_low+1); %coarse binning
bins_az_low = linspace(prctile(RF_parameters_2D(:,4),1),prctile(RF_parameters_2D(:,4),99),nbins_low+1);%coarse binning
bins_el_high = linspace(elmin,elmax,nbins_high+1);%fine binning
bins_az_high = linspace(prctile(RF_parameters_2D(:,4),1),prctile(RF_parameters_2D(:,4),99),nbins_high+1);%fine binning
% elevation and azimuth bin centers
bincenters_el_low=bins_el_low(1:end-1)+median(diff(bins_el_low))/2;
bincenters_az_low=bins_az_low(1:end-1)+median(diff(bins_az_low))/2;
bincenters_el_high=bins_el_high(1:end-1)+median(diff(bins_el_high))/2;
bincenters_az_high=bins_az_high(1:end-1)+median(diff(bins_az_high))/2;

% make 1-d-RF binned averages with coarse binning
[N_AzEl1d,N_El1d,N_Az1d,~,avg1dRF_el,avg1dRF_az,~,avg1dRF_params_el_az,~,avg1dRF_params_el,~,avg1dRF_params_az]=fit_mean_RFs(RF_1D,RF_parameters_2D,bins_az_low,bins_el_low,bins_az_low,bins_el_low,1);
% extract relevant parameters from 2D binning (center size, surround strength, surround assymetry) 
avg1dRF_params_el_az=cell2mat(cellfun(@(x) shiftdim(x,-1),avg1dRF_params_el_az(:,:,end),'UniformOutput',false));
avg1dRF_params_el_az=avg1dRF_params_el_az(:,:,[2 1 3]);
SEMAz=std(avg1dRF_params_el_az,0,1,'omitnan')./sqrt(sum(~isnan(avg1dRF_params_el_az),1));
SEMEl=std(avg1dRF_params_el_az,0,2,'omitnan')./sqrt(sum(~isnan(avg1dRF_params_el_az),2));
% extract relevant parameters from 1D elevation binning (center size, surround strength, surround assymetry)
avg1dRF_params_el=cell2mat(cellfun(@(x) shiftdim(x,-1),avg1dRF_params_el(:,:,end),'UniformOutput',false));
avg1dRF_params_el=avg1dRF_params_el(:,:,[2 1 3]);
LowSampling=repmat(N_El1d(:,:,end)<minRGCparam,1,1,size(avg1dRF_params_el,3));
avg1dRF_params_el(LowSampling)=nan;
% extract relevant parameters from 1D azimuth binning (center size, surround strength, surround assymetry)
avg1dRF_params_az=cell2mat(cellfun(@(x) shiftdim(x,-1),avg1dRF_params_az(:,:,end),'UniformOutput',false));
avg1dRF_params_az=avg1dRF_params_az(:,:,[2 1 3]);
LowSampling=repmat(N_Az1d(:,:,end)<minRGCparam,1,1,size(avg1dRF_params_az,3));
avg1dRF_params_az(LowSampling)=nan;

% make 1-d-RF binned averages with fine binning
[~,~,~,~,avg1dRF_fine_el,avg1dRF_fine_az]=fit_mean_RFs(RF_1D,RF_parameters_2D,bins_az_high,bins_el_high,[],[],0);

% make 2-d-RF binned averages
[N_AzEl2d,N_El2d,N_Az2d,avg2dRF_el_az,avg2dRF_el,avg2dRF_az]=fit_mean_RFs(RF_2D,RF_parameters_2D,bins_az_low,bins_el_low,bins_az_low,bins_el_low,0);


%% %%%%%%%%%%%%%%%%%%%%%%%%% STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ks test of azimuthal bin values above & below horizon
stat_p_el=nan(size(avg1dRF_params_el_az,3),1);
for sb=1:size(avg1dRF_params_el_az,3)
    sbb=discretize(bincenters_el_low,statbins_el);
    temp1=reshape(avg1dRF_params_el_az(sbb==1,:,sb),[],1);
    temp2=reshape(avg1dRF_params_el_az(sbb==3,:,sb),[],1);
    [~,stat_p_el(sb)]=kstest2(temp1,temp2);
    
end
% ks test of elevation bin values nasal & temporal of 65°
stat_p_az=nan(size(avg1dRF_params_el_az,3),1);
for sb=1:size(avg1dRF_params_el_az,3)
    sbb=discretize(bincenters_az_low,statbins_az);
    temp1=reshape(avg1dRF_params_el_az(:,sbb==1,sb),[],1);
    temp2=reshape(avg1dRF_params_el_az(:,sbb==3,sb),[],1);
    [~,stat_p_az(sb)]=kstest2(temp1,temp2);
    
end
% Linear regression
lr=nan(2,3,2);
ss=RF_parameters_2D(:,3)>=bins_el_low(1) & RF_parameters_2D(:,3)<=bins_el_low(end) & RF_parameters_2D(:,4)>=bins_az_low(1) & RF_parameters_2D(:,4)<=bins_az_low(end);
for x=1:3
    mdl=fitlm(RF_parameters_2D(ss,3:4),RF_parameters_1D(ss,x),'VarNames',{'Elevation','Azimuth','param'});
    lr(1,x,1)=mdl.Coefficients{'Elevation','Estimate'};
    lr(1,x,2)=mdl.Coefficients{'Azimuth','Estimate'};
    lr(2,x,1)=mdl.Coefficients{'Elevation','pValue'};
    lr(2,x,2)=mdl.Coefficients{'Azimuth','pValue'};
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 5 d (Example RF_2D elevation)
crange=[-.5 .5];
selx=round(linspace(1,numel(bincenters_el_low),min(numel(bincenters_el_low),8))); % sample across range
fig_exp_el=figure('Position',[ 55          60        1824         993],'Color','white') ;
tt=tiledlayout(numel(selx),16);
title(tt,'Example RF_2D','Elevation')
for nx=numel(selx):-1:1
    sx=find(abs(RF_parameters_2D(:,3)-bincenters_el_low(selx(nx)))<median(diff(bincenters_el_low))/4);
    [~,sa]=sort(RF_SNR_1D(sx),'descend');
    sxx=sx(sa);
    for nxx=1:2:16
        
        nexttile;
        if nxx<numel(sxx)
        imagesc(RF_2D(:,:,sxx(nxx)),crange); colormap(flipud(othercolor('RdBu11')));
        axis xy equal;
        xlim([0 50]);ylim([0 50])
        set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
        title(sprintf('el: %.1f deg',RF_parameters_2D(sxx(nxx),3)))
        end
        nexttile;
        
        if nxx<numel(sxx)
        plot(RF_1D(:,:,sxx(nxx)),-25:25,'k')
        xlim([-1 1])
        ylim([-25 25])
        set(gca,'Color','none');
        box off
        ylabel('relative RF el (°)')
        end
    end
end
exportgraphics(fig_exp_el,[figdir 'exampleRF_elevation.pdf'],'BackgroundColor','none','ContentType','vector')

%% Fig.5e & supplementary Fig. 9h (1D RF visualization)

% filter width in pixels
saz=round(mean_filt_width/median(diff(bincenters_az_high)));
sel=round(mean_filt_width/median(diff(bincenters_el_high)));

fig_1d_vis=figure('Position',[21         414        1851         657]);
tt=tiledlayout(2,1);
title(tt,'1D meanRF visualization')

nt(1)=nexttile(1);
bin_indicators=bincenters_el_low'+[-median(diff(bincenters_el_low)) median(diff(bincenters_el_low))]/2;
RFelAvg1n=avg1dRF_fine_el./max(abs(avg1dRF_fine_el),[],1);
imagesc(bincenters_el_high,-25:25,smoothdata(permute(RFelAvg1n(:,1,:,1,end),[1 3 2]),2,'movmean',sel),crange);axis xy; colormap(flipud(othercolor('RdBu11')));
ylim([-18 18])
xl=xlim;
set(gca,'XTick',[xl(1) 0 xl(end)]);
hold on;
line(bin_indicators(1:2:end-1,:)',repmat(-17,size(bin_indicators(1:2:end-1,:)))','Color','k')
line(bin_indicators(2:2:end,:)',repmat(-16,size(bin_indicators(2:2:end,:)))','Color','k')

title('Elevation',sprintf('binsize=%.2f',median(diff(bincenters_el_high))))
ylabel('1D relative RF elevation (deg)')
xlabel('RF center elevation (deg)')
box off

nt(2)=nexttile(2);
bin_indicators=bincenters_az_low'+[-median(diff(bincenters_az_low)) median(diff(bincenters_az_low))]/2;
RFazAvg1n=avg1dRF_fine_az./max(abs(avg1dRF_fine_az),[],1);
imagesc(bincenters_az_high,-25:25,smoothdata(permute(RFazAvg1n(:,1,1,:,end),[1 4 3 2]),2,'movmean',saz),crange);axis xy; colormap(flipud(othercolor('RdBu11')));
ylim([-18 18])
hold on;
line(bin_indicators(1:2:end-1,:)',repmat(-17,size(bin_indicators(1:2:end-1,:)))','Color','k')
line(bin_indicators(2:2:end,:)',repmat(-16,size(bin_indicators(2:2:end,:)))','Color','k')
xl=xlim;
set(gca,'XTick',[xl(1) xl(end)]);
title('Azimuth',sprintf('binsize=%.2f',median(diff(bincenters_az_high))))
ylabel('1D relative RF elevation (deg)')
xlabel('RF center azimuth (deg)')
box off
exportgraphics(fig_1d_vis,[figdir 'tunnel_plots_el_az.pdf'],'BackgroundColor','none','ContentType','vector')

%% Fig 5f (mean RF_2D elevation)
crange=[-.3 .3];
selx=round(linspace(1,numel(bincenters_el_low),min(numel(bincenters_el_low),inf)));
fig_mean_binned_el=figure('Position',[ 55          61        1704         992],'Color','white') ;
tt=tiledlayout(numel(selx)/4,8);
title(tt,'Average RF_2D','Elevation')
for nx=1:numel(selx)
    nexttile;
    imagesc(avg2dRF_el(:,:,selx(nx),1,end),crange); colormap(flipud(othercolor('RdBu11')));
    axis xy equal;
    xlim([0 50]);ylim([0 50])
    set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
    title(sprintf('el: %.1f deg, N=%u',bincenters_el_low(selx(nx)),N_El2d(selx(nx),1,end)))
    nexttile;
    plot(avg1dRF_el(:,:,selx(nx),1,end),-25:25,'k')
    xlim([-1.1 .5])
    ylim([-25 25])
    set(gca,'Color','none');
    box off
    ylabel('relative RF el (°)')
end
exportgraphics(fig_mean_binned_el,[figdir 'avgRF_elevation.pdf'],'BackgroundColor','none','ContentType','vector')
%

%% Fig.5 g,h,i (1D elevation parameters)

fig_mean_binned_el1D=figure('Position',[21         726        1715         345]);
tt=tiledlayout(1,15);
title(tt,'1D meanRF parameterization','Elevation')
fig_mean_binned_el1D=plot_1D_params(repmat(statbins_el,3,1),stat_p_el,avg1dRF_params_el,SEMEl,cat(2,bincenters_el_low,fliplr(bincenters_el_low)),bincenters_el_low,'RF center elevation (deg)',lr,fig_mean_binned_el1D);
exportgraphics(fig_mean_binned_el1D,[figdir '1D_elevation_parameterization_1D_avgRFs.pdf'],'BackgroundColor','none','ContentType','vector')

%% supplementary Fig 9 b,c (saccade axis example and all recordings)
example_num=5;
saccade_axis_angle_rad=deg2rad(saccade_axis(:,1));

fig_saccade_axis_all=figure('Position',[416   667   825   432]);
tiledlayout(1,2);
nexttile;
polarscatter(saccade_parameters{example_num}(:,2),saccade_parameters{example_num}(:,3),20,'k.');hold on
polarplot([saccade_axis_angle_rad(example_num) saccade_axis_angle_rad(example_num)+pi]',max(saccade_parameters{example_num}(:,3))*ones(numel(saccade_axis_angle_rad(example_num)),2)','r'); hold off
title('Saccade axes',sprintf('saccade axis tuning= %.2f',mean(saccade_axis(example_num,2))))

nexttile;
polarplot([saccade_axis_angle_rad saccade_axis_angle_rad+pi]',1*ones(numel(saccade_axis_angle_rad),2)','k'); hold on
polarplot([saccade_axis_angle_rad(example_num) saccade_axis_angle_rad(example_num)+pi]',1*ones(numel(saccade_axis_angle_rad(example_num)),2)','r'); hold off
rlim([0 1])
title('Saccade axes',sprintf('mean saccade axis tuning= %.2f',mean(saccade_axis(:,2))))
exportgraphics(fig_saccade_axis_all,[figdir 'saccade_axes_all.pdf'],'BackgroundColor','none','ContentType','vector')


%% supplementary Fig 9 d,e,f (2D parameterization)

AD=double(N_AzEl1d(:,:,end)>=5);
fig_binned_val_1davg=figure('Position',[32         470        1654         549],'Color','white') ;
tt=tiledlayout(1,3);
title(tt,'2D parameterization of 1D Average RF_2D')
nt(1)=nexttile(1);
imagesc(bincenters_az_low,bincenters_el_low,avg1dRF_params_el_az(:,:,1),'AlphaData',AD,[0 .6]);axis xy equal
colormap(nt(1),flipud(othercolor('PRGn11')))
nt(1).Color=[1 1 1];
cb(1)=colorbar;
xl=xlim;
yl=ylim;
set(gca,'XTick',[xl(1) xl(end)],'YTick',[yl(1) 0 yl(end)]);

ylabel('el vis angle °')
xlabel('az vis angle °')
title('relative surround strength')

nt(2)=nexttile(2);
imagesc(bincenters_az_low,bincenters_el_low,avg1dRF_params_el_az(:,:,2),'AlphaData',AD,[7 10]);axis xy equal
colormap(nt(2),flipud(othercolor('PRGn11')))
nt(2).Color=[1 1 1];
cb(2)=colorbar;
ylabel('el vis angle °')
xlabel('az vis angle °')
title('center size')
xl=xlim;
yl=ylim;
set(gca,'XTick',[xl(1) xl(end)],'YTick',[yl(1) 0 yl(end)]);

nt(3)=nexttile(3);
imagesc(bincenters_az_low,bincenters_el_low,avg1dRF_params_el_az(:,:,3),'AlphaData',AD,[-1 1]);axis xy equal
colormap(nt(3),othercolor('BrBG11'))
nt(3).Color=[1 1 1];
cb(3)=colorbar;
ylabel('el vis angle °')
xlabel('az vis angle °')
title('vertical surround asymmetry')
xl=xlim;
yl=ylim;
set(gca,'XTick',[xl(1) xl(end)],'YTick',[yl(1) 0 yl(end)]);

exportgraphics(fig_binned_val_1davg,[figdir '2D_parameterization_1D_avgRFs.pdf'],'BackgroundColor','none','ContentType','vector')
%% supplementary Fig 9g (example RF_2D azimuth)
crange=[-.5 .5];
selx=round(linspace(1,numel(bincenters_az_low),min(numel(bincenters_az_low),8)));
fig_exp_az=figure('Position',[ 55          60        1824         993],'Color','white') ;
tt=tiledlayout(numel(selx),16);
title(tt,'Example RF_2D','Azimuth')
for nx=numel(selx):-1:1
    sx=find(abs(RF_parameters_2D(:,4)-bincenters_az_low(selx(nx)))<median(diff(bincenters_az_low))/4);
    [~,sa]=sort(RF_SNR_1D(sx),'descend');
    sxx=sx(sa);
    for nxx=1:2:16
        
        nexttile;
        imagesc(RF_2D(:,:,sxx(nxx)),crange); colormap(flipud(othercolor('RdBu11')));
        axis xy equal;
        xlim([0 50]);ylim([0 50])
        set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
        title(sprintf('az: %.1f deg',RF_parameters_2D(sxx(nxx),4)))
        
        nexttile;
        plot(RF_1D(:,:,sxx(nxx)),-25:25,'k')
        xlim([-1 1])
        ylim([-25 25])
        set(gca,'Color','none');
        box off
        ylabel('relative RF el (°)')
    end
end
exportgraphics(fig_exp_az,[figdir 'exampleRF_azimuth.pdf'],'BackgroundColor','none','ContentType','vector')

%% supplementary Fig 9i (mean RF azimuth)
selx=round(linspace(1,numel(bincenters_az_low),min(numel(bincenters_az_low),inf)));
fig_mean_binned_az=figure('Position',[ 55          61        1704         992],'Color','white') ;
tt=tiledlayout(numel(selx)/4,8);
title(tt,sprintf('Average RF_2D, binsize=%.1f°',median(diff(bincenters_az_low))),'Azimuth')
for nx= 1:numel(selx)
    nexttile;
    imagesc(avg2dRF_az(:,:,1,selx(nx),end),crange); colormap(flipud(othercolor('RdBu11')));
    axis xy equal;
    box off
    set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
    xlim([0 50]);ylim([0 50])
    title(sprintf('az: %.1f deg, N=%u',bincenters_az_low(selx(nx)),N_Az2d(1,selx(nx),end)))
    nexttile;%(nx+numel(selx));
    plot(avg1dRF_az(:,:,1,selx(nx),end),-25:25,'k')
    set(gca,'Color','none');
    ylabel('relative RF el (°)')
    box off
    xlim([-1.1 .5])
    ylim([-25 25])
end
exportgraphics(fig_mean_binned_az,[figdir 'avgRF_azimuth.pdf'],'BackgroundColor','none','ContentType','vector')
%% supplementary Fig.10 j,k,l (1d parameters azimuth)

curVal=permute(avg1dRF_params_az,[2 1 3]);
curSEM=permute(SEMAz,[2 1 3]);
cur_xsem=cat(2,bincenters_az_low,fliplr(bincenters_az_low));
cur_x=bincenters_az_low;
cur_xlabel='RF center azimuth (deg)';
fig_mean_binned_az1D=figure('Position',[22         292        1715         345]);
tt=tiledlayout(1,15);
title(tt,sprintf('1D meanRF parameterization, binsize=%.1f°',median(diff(bincenters_az_low))),'Azimuth')
fig_mean_binned_az1D=plot_1D_params(repmat(statbins_az,3,1),stat_p_az,curVal,curSEM, cur_xsem, cur_x,cur_xlabel, lr,fig_mean_binned_az1D);

exportgraphics(fig_mean_binned_az1D,[figdir '1D_azimuth_parameterization_1D_avgRFs.pdf'],'BackgroundColor','none','ContentType','vector')


%% supplementary Fig 5b (2D avg RF)
crange=[-.3 .3];

fig_mean_binned_both=figure('Position',[50          50        1056         852],'Color','white') ;

tt=tiledlayout(numel(bincenters_el_low),numel(bincenters_az_low),'TileSpacing','none','Padding','none');
nxx=0;
title(tt,sprintf('Average RF_2D, binsEl=%.1f:%.1f:%.1f°, binsAz=%.1f:%.1f:%.1f°',bincenters_el_low(1),median(diff(bincenters_el_low)),bincenters_el_low(end),bincenters_az_low(1),median(diff(bincenters_az_low)),bincenters_az_low(end)),'Elevation & Azimuth')
for nxe=numel(bincenters_el_low):-1:1
    for nxa=1:numel(bincenters_az_low)
        nxx=nxx+1;
        nexttile(nxx);
        if N_AzEl2d(nxe,nxa,end)>=20
            imagesc(avg2dRF_el_az(:,:,nxe,nxa,end),crange);axis xy equal; colormap(flipud(othercolor('RdBu11')));
            box off
            set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
        end
        box off
        set(gca,'XTick',[],'YTick',[],'Color','none','YColor','none','XColor','none');
    end
end
exportgraphics(fig_mean_binned_both,[figdir 'avgRF_el_az.pdf'],'BackgroundColor','none','ContentType','vector')


%% %%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [N_AzEl,N_El,N_Az,RF2avg,RFelAvg,RFazAvg,RFprops_AzElAvg,RFfparams_AzElAvg,RFprops_ElAvg,RFfparams_ElAvg,RFprops_AzAvg,RFfparams_AzAvg]=fit_mean_RFs(RF,RF_parameters_2D,binsAz1,binsEl1,binsAz2,binsEl2,ftd)

%pre allocate
RFprops_AzElAvg=[];
RFfparams_AzElAvg=[];
RFprops_ElAvg=[];
RFfparams_ElAvg=[];
RFprops_AzAvg=[];
RFfparams_AzAvg=[];
RF2avg=nan([size(RF,[1 2]) numel(binsEl2)-1 numel(binsAz2)-1]);
RFelAvg=nan([size(RF,[1 2])  numel(binsEl1)-1 1]);
RFazAvg=nan([size(RF,[1 2])  1 numel(binsAz1)-1]);
N_AzEl=nan([numel(binsEl2)-1 numel(binsAz2)-1]);
N_El=nan([numel(binsEl1)-1 1]);
N_Az=nan([1 numel(binsAz1)-1]);

% bin locations
if ~isempty(binsEl1);El_idx1 = discretize(RF_parameters_2D(:,3),binsEl1);end
if ~isempty(binsAz1);Az_idx1 = discretize(RF_parameters_2D(:,4),binsAz1);end
if ~isempty(binsEl2);El_idx2 = discretize(RF_parameters_2D(:,3),binsEl2);end
if ~isempty(binsAz2);Az_idx2 = discretize(RF_parameters_2D(:,4),binsAz2);end


% compute mean RF_2D over 2d space
if ~isempty(binsEl2) && ~isempty(binsAz2)
    for nxe=1:numel(binsEl2)-1
        for nxa=1:numel(binsAz2)-1

            RF2avg(:,:,nxe,nxa)=mean(RF(:,:,El_idx2==nxe & Az_idx2==nxa),3); %mean RF
            N_AzEl(nxe,nxa)=sum(El_idx2==nxe & Az_idx2==nxa); % number of boutons
        end
    end
    RF2avg=RF2avg./max(abs(RF2avg),[],[1 2]); %normalize
    [RFprops_AzElAvg,RFfparams_AzElAvg]=fit_avg_RF(RF2avg,ftd); %fit
end

% compute mean RF_2D over elevation
if ~isempty(binsEl1)
    for nxe=1:numel(binsEl1)-1
        
        RFelAvg(:,:,nxe,1)=mean(RF(:,:,El_idx1==nxe),3); %mean RF
        N_El(nxe,1)=sum(El_idx1==nxe);% number of boutons
    end
    RFelAvg=RFelAvg./max(abs(RFelAvg),[],[1 2]); %normalize
    [RFprops_ElAvg,RFfparams_ElAvg]=fit_avg_RF(RFelAvg,ftd);%fit
end

% compute mean RF_2D over azimuth
if ~isempty(binsAz1)
    for nxa=1:numel(binsAz1)-1

        RFazAvg(:,:,1,nxa)=mean(RF(:,:,Az_idx1==nxa),3); %mean RF
        N_Az(1,nxa)=sum(Az_idx1==nxa);% number of boutons
    end
    RFazAvg=RFazAvg./max(abs(RFazAvg),[],[1 2]); %normalize
    [RFprops_AzAvg,RFfparams_AzAvg]=fit_avg_RF(RFazAvg,ftd);%fit
end

end

%%
function [RFprops,RFfparams]=fit_avg_RF(RFavg,RF_dimensions)

RFprops=cell(size(RFavg,3),size(RFavg,4),size(RFavg,5));
RFfparams=cell(size(RFavg,3),size(RFavg,4),size(RFavg,5));

if RF_dimensions>0
    NN=numel(RFprops(:,:,1));nn=0;
    w=waitbar(0,'fitting mean RF_2D');
    for nx1=1:size(RFprops,1)
        for nx2=1:size(RFprops,2)
            nn=nn+1;
            tempp=double(RFavg(:,:,nx1,nx2,:));
            for nx3=size(RFprops,3)
                if ~all(isnan(tempp(:,:,:,:,nx3)))
                    if RF_dimensions==1
                        test=tempp(:,:,:,:,nx3);
                        [pH,pL,pW]=findpeaks(-test,'NPeaks',1,'MinPeakHeight',max(-test)*.9); %init peak
                        % set init values and bounds
                        if ~isempty(pH)
                            pW=2*sqrt(pW);
                            StartPoint=[pH,     pL,     pW,     0,      pL,             10*pW];
                            Lower=     [pH/2,   pL-20,  .25*pW, 0,      0,              pW];
                            Upper=     [inf,    pL+20,  inf,    2*pH,   numel(test),    inf];
                        else
                            [pH,pL]=max(-test);
                            StartPoint=[pH,     pL,     5,     0,   pL,    10];
                            Lower=     [pH/2,     0,  1, 0,   0,     1];
                            Upper=     [inf,   numel(test),  inf,   pH,  numel(test),   inf];
                            
                        end
                        RF_fittype_1d =  fittype('amp1*exp(-((x-x01).^2/(2*wx1^2)))-amp2*exp(-((x-x02).^2/(2*wx2^2)))', 'dependent',{'z'},'independent',{'x'},...
                            'coefficients',{'amp1','x01','wx1','amp2','x02','wx2'});
                        RF_fitoptions=fitoptions(RF_fittype_1d);
                        RF_fitoptions.StartPoint=StartPoint;
                        RF_fitoptions.Lower= Lower;
                        RF_fitoptions.Upper= Upper;
                        [RF_fitresult, RF_gof] = fit( (1:size(test,1))', -test, RF_fittype_1d, RF_fitoptions);
                        elvals=(1:size(test,1))';
                        elvals2=linspace(1,size(test,1),10*numel(elvals))';
                        model=RF_fitresult(elvals);
                        centvals=test.*(model>0);
                        survals=test.*(model<0);
                        com=[sum(elvals.*centvals)/sum(centvals) sum(elvals.*survals)/sum(survals)];
                        half_size=floor(size(test,1)/2);
                        center_mask=false(size(test));
                        center_mask(elvals>=RF_fitresult.x01-2*RF_fitresult.wx1 & elvals<=RF_fitresult.x01+2*RF_fitresult.wx1)=true;
                        center_mask2=false(size(elvals2));
                        center_mask2(elvals2>=RF_fitresult.x01-2*RF_fitresult.wx1 & elvals2<=RF_fitresult.x01+2*RF_fitresult.wx1)=true;
                        surround_mask=false(size(test));
                        surround_mask(elvals>=RF_fitresult.x02-2*RF_fitresult.wx2 & elvals<=RF_fitresult.x02+2*RF_fitresult.wx2)=true;
                        surround_mask(center_mask)=false;
                        surround_mask_upper=surround_mask;
                        surround_mask_upper2=surround_mask;
                        surround_mask_upper(1:half_size)=false;
                        surround_mask_upper2(1:floor(RF_fitresult.x01))=false;
                        surround_mask_lower=surround_mask;
                        surround_mask_lower2=surround_mask;
                        surround_mask_lower(end-half_size:end)=false;
                        surround_mask_lower2(ceil(RF_fitresult.x01):end)=false;
                        sur_sum=abs(sum(test(surround_mask)));
                        cent_sum=abs(sum(test(center_mask)));
                        rel_surround_strength=abs(sur_sum./cent_sum);
                        sur_sum_l=abs(sum(test(surround_mask_lower)));
                        sur_sum_u=abs(sum(test(surround_mask_upper)));
                        sur_sum_l2=abs(sum(test(surround_mask_lower2)));
                        sur_sum_u2=abs(sum(test(surround_mask_upper2)));
                        vertical_surround_assymetry=(sur_sum_u-sur_sum_l)./(sur_sum_u+sur_sum_l);
                        vertical_surround_assymetry2=(sur_sum_u2-sur_sum_l2)./(sur_sum_u2+sur_sum_l2);
                        center_size=2*RF_fitresult.wx1;
                        center_surround_distance=abs(diff(com));
                        RFprops{nx1,nx2,nx3}=RF_fitresult;
                        RFfparams{nx1,nx2,nx3}=cat(2,center_size,rel_surround_strength,vertical_surround_assymetry,center_surround_distance,RF_gof.rsquare,sum(center_mask2)/10,vertical_surround_assymetry2);
                    else
                        test=tempp(:,:,:,:,nx3);
                        RFprops{nx1,nx2,nx3} = RF_gaussian_mixture_spatial_only(test);
                        R2vals=mk_R2_RF(RFprops{nx1,nx2,nx3});
                        
                        center_mask=logical( RFprops{nx1,nx2,nx3}.center_mask);
                        surround_mask=logical( RFprops{nx1,nx2,nx3}.surround_mask);
                        half_size=size(test,1)/2;
                        surround_mask_upper=surround_mask;
                        surround_mask_upper(1:floor(half_size),:)=false;
                        surround_mask_lower=surround_mask;
                        surround_mask_lower(ceil(half_size):end,:)=false;
                        
                        surRF=test;
                        surRF(~surround_mask)=0;
                        surRF(center_mask)=0;
                        sur_sum=abs(squeeze(sum(surRF,[1 2])));
                        centRF=test;
                        centRF(~center_mask)=0;
                        cent_sum=abs(squeeze(sum(centRF,[1 2])));
                        rel_surround_strength=abs(sur_sum./cent_sum);
                        
                        sur_sum_l=abs(sum(test(surround_mask_lower)));
                        sur_sum_u=abs(sum(test(surround_mask_upper)));
                        
                        vertical_surround_assymetry=(sur_sum_u-sur_sum_l)./(sur_sum_u+sur_sum_l);
                        vertical_surround_assymetry2=nan;
                        center_size=2*sqrt(sum(center_mask,'all')./(2*pi));
                        center_size2=2*sqrt(sum(center_mask,'all')./(2*pi));
                        center_surround_distance=RFprops{nx1,nx2,nx3}.distance_btw_centers_of_mass;
                        RFfparams{nx1,nx2,nx3}=cat(2,center_size,rel_surround_strength,vertical_surround_assymetry,center_surround_distance,center_size2,R2vals.R2_cs,vertical_surround_assymetry2);
                    end
                else
                    RFfparams{nx1,nx2,nx3}=nan(1,7);
                end
            end
            waitbar(nn/NN,w)
        end
    end
    close(w)
end
end

%%
function fig=plot_1D_params(statbins,stat_p,curVal,curSEM,cur_xsem,cur_x,cur_xlabel,lr,fig)

m=mean(statbins(:,2:3),2);
bcc=categorical({'Elevation','Azimuth'},{'Elevation','Azimuth'});
nt(1)=nexttile(1,[1 3]);
sem_plot=cat(1,curVal(:,1,1)+curSEM(:,1,1),flipud(curVal(:,1,1)-curSEM(:,1,1)))';
fill(cur_xsem(~isnan(sem_plot)), sem_plot(~isnan(sem_plot)),'k','EdgeColor','none','FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on
plot(cur_x,curVal(:,1,1),'Color',[0 0 0]);
xlim(cur_x([1 end]));
if cur_x(1) < 0
    set(gca,'XTick',[cur_x(1) 0 cur_x(end)]);
else
    set(gca,'XTick',[cur_x(1) cur_x(end)]);
end
xlabel(cur_xlabel)
title('relative surround strength')
ylim([0 0.8]);
line([min(cur_x) statbins(1,2);statbins(1,3) max(cur_x)]',repmat(0.75,2,2),'Color','k')
text(m(1),0.75,sprintf('P=%.2e',stat_p(1)),'HorizontalAlignment','center','VerticalAlignment','bottom')
box off

bt(1)=nexttile(5);
b=bar(bcc,[lr(1,1,1) lr(1,1,2)]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1(1),ytips1(1),sprintf('p=%.1e',lr(2,1,1)),'HorizontalAlignment','center','VerticalAlignment','bottom')
text(xtips1(2),ytips1(2),sprintf('p=%.1e',lr(2,1,2)),'HorizontalAlignment','center','VerticalAlignment','bottom')

title('Lin.Reg. Coefs')
box off

nt(2)=nexttile(6,[1 3]);
sem_plot=cat(1,curVal(:,1,2)+curSEM(:,1,2),flipud(curVal(:,1,2)-curSEM(:,1,2)))';
fill(cur_xsem(~isnan(sem_plot)), sem_plot(~isnan(sem_plot)),'k','EdgeColor','none','FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on
plot(cur_x,curVal(:,1,2),'Color',[0 0 0]);
xlim(cur_x([1 end]));
xlabel(cur_xlabel)
if cur_x(1) < 0
    set(gca,'XTick',[cur_x(1) 0 cur_x(end)]);
else
    set(gca,'XTick',[cur_x(1) cur_x(end)]);
end
ylabel('diameter (deg)')
title('center size')
line([min(cur_x) statbins(2,2);statbins(2,3) max(cur_x)]',repmat(9.75,2,2),'Color','k')
text(m(2),9.75,sprintf('P=%.2e',stat_p(2)),'HorizontalAlignment','center','VerticalAlignment','bottom')
box off


bt(2)=nexttile(10);
b=bar(bcc,[lr(1,2,1) lr(1,2,2)]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1(1),ytips1(1),sprintf('p=%.1e',lr(2,2,1)),'HorizontalAlignment','center','VerticalAlignment','bottom')
text(xtips1(2),ytips1(2),sprintf('p=%.1e',lr(2,2,2)),'HorizontalAlignment','center','VerticalAlignment','bottom')

title('Lin.Reg. Coefs')
box off

nt(3)=nexttile(11,[1 3]);
sem_plot=cat(1,curVal(:,1,3)+curSEM(:,1,3),flipud(curVal(:,1,3)-curSEM(:,1,3)))';
fill(cur_xsem(~isnan(sem_plot)), sem_plot(~isnan(sem_plot)),'k','EdgeColor','none','FaceColor',[.2 .2 .2],'FaceAlpha',.5);hold on
plot(cur_x,curVal(:,1,3),'Color',[0 0 0]);
xlim(cur_x([1 end]));
if cur_x(1) < 0
    set(gca,'XTick',[cur_x(1) 0 cur_x(end)]);
else
    set(gca,'XTick',[cur_x(1) cur_x(end)]);
end
ylim([-.25 1.25]);
xlabel(cur_xlabel)
title('vertical surround asymmetry')
line([min(cur_x) statbins(3,2);statbins(3,3) max(cur_x)]',repmat(1.2,2,2),'Color','k')
text(m(3),1.2,sprintf('P=%.2e',stat_p(3)),'HorizontalAlignment','center','VerticalAlignment','bottom')
box off

bt(3)=nexttile(15);
b=bar(bcc,[lr(1,3,1) lr(1,3,2)]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
text(xtips1(1),ytips1(1),sprintf('p=%.1e',lr(2,3,1)),'HorizontalAlignment','center','VerticalAlignment','bottom')
text(xtips1(2),ytips1(2),sprintf('p=%.1e',lr(2,3,2)),'HorizontalAlignment','center','VerticalAlignment','bottom')

title('Lin.Reg. Coefs')
box off

end
