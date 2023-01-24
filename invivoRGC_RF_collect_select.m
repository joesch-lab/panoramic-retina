
%% script to collect receptive field estimates and parameterization
resultfiles=dir('%path to files%\in_vivo_RF_data\*Result.mat');

collected_savename='D:\RGC_invivo_RFs_collected.mat';
selected_savename='D:\RGC_invivo_RFs_selected.mat';

%% collect data across recordings
RFa_bgA={};
RFProps_a={};
RFa_param=cell(numel(resultfiles),1);
sac_deg_tuningA=nan(numel(resultfiles),2);
saccade_params_a=cell(numel(resultfiles),1);
%
for rf=1:numel(resultfiles)

    clear('saccade_params','sac_deg_tuning','RF_param','paramsA','GeneralParams','RFa_bg');
    %load
    rfile=fullfile(resultfiles(rf).folder,resultfiles(rf).name);
    fprintf('load file %u\n',rf)
    load(rfile,'RF_param_key','RFprops','saccade_params','sac_deg_tuning','RF_param','paramsA','GeneralParams','RFa_bg');
    sac_deg_tuningA(rf,:)=sac_deg_tuning;
    saccade_params_a{rf}=saccade_params;
    RFa_bg_norot=RFa_bg;
    
    rotA=sac_deg_tuningA(rf,1);
    %rotate RF center positions according to saccade axis
    if ~isnan(rotA)
        R = makehgtform('yrotate',deg2rad(-rotA)); %rotation matrix
        R = R(1:3,1:3);
        coords=nan(size(RF_param,1),3);
        [coords(:,1),coords(:,2),coords(:,3)] = sph2cart(deg2rad(RF_param(:,25)+25),deg2rad(RF_param(:,24)),40) ;%spherical=>cartesian coordinates with eye axis
        rotcoords_cart=coords;
        for c=1:size(coords,1)
            rotcoords_cart(c,:) = (R*coords(c,:)')';%rotate
        end
        rotcoords_sph=rotcoords_cart;
        [rotcoords_sph(:,1),rotcoords_sph(:,2),rotcoords_sph(:,3)]=cart2sph(rotcoords_cart(:,1),rotcoords_cart(:,2),rotcoords_cart(:,3));%back to spherical coordinates
        RF_param(:,3:4)=rad2deg(rotcoords_sph(:,[2 1]))-[0 25];
        
    end
    
    w=waitbar(0);
    for n=1:size(RF_param,1)
        %extract parameters
        RFs_bg=RFa_bg{n};
        [~,RF_param(n,8)]=max(std(RFs_bg,0,[1 2]));% max deviation
        RF_param(n,9)=rf; % recording number
        RF_param(n,10)=n; % roi number
        
        if ~isnan(rotA)
            temp=rot90(RFs_bg,2);
            
            for t=1:size(temp,4)
                temp(:,:,1,t)=imrotate(temp(:,:,1,t),-rotA,'bilinear','crop'); %rotate RF
                
            end
            RFs_bg=rot90(temp,2);
        end
        RFa_bg{n}=RFs_bg;
        
        
        if RFprops{n}.center_exist==1
            %rotate masks
            if ~isnan(rotA)
                if ~isfield(RFprops{n},'center_mask_norot')
                    RFprops{n}.center_mask_norot=RFprops{n}.center_mask;
                    RFprops{n}.surround_mask_norot=RFprops{n}.surround_mask;
                else
                    RFprops{n}.center_mask=RFprops{n}.center_mask_norot;
                    RFprops{n}.surround_mask=RFprops{n}.surround_mask_norot;
                end
                RFprops{n}.center_mask=rot90(imrotate(rot90( RFprops{n}.center_mask,2),-rotA,'nearest','crop'),2);
                RFprops{n}.surround_mask=rot90(imrotate(rot90( RFprops{n}.surround_mask,2),-rotA,'nearest','crop'),2);
            end
            %compute RF parameterizations
            center_mask=logical(RFprops{n}.center_mask);
            surround_mask=logical(RFprops{n}.surround_mask);
            surround_mask(center_mask)=false;
            half_size=size(RFs_bg,1)/2;
            surround_mask_upper=surround_mask;
            surround_mask_upper(1:floor(half_size),:)=false;
            surround_mask_lower=surround_mask;
            surround_mask_lower(ceil(half_size):end,:)=false;
            surRF=RFs_bg(:,:,:,RF_param(n,8));
            surRF(~surround_mask)=0;
            sur_sum=abs(sum(surRF,'all'));
            centRF=RFs_bg(:,:,:,RF_param(n,8));
            centRF(~center_mask)=0;
            cent_sum=abs(sum(centRF,'all'));
            rel_surround_strength=abs(sur_sum./cent_sum);
            temp=RFs_bg(:,:,:,RF_param(n,8));
            sur_sum_l=abs(sum(temp(surround_mask_lower)));
            sur_sum_u=abs(sum(temp(surround_mask_upper)));
            vertical_surround_assymetry=(sur_sum_u-sur_sum_l)./(sur_sum_u+sur_sum_l);
            sur_siz_l=sum(surround_mask_lower,'all');
            sur_siz_u=sum(surround_mask_upper,'all');
            vertical_surround_size_assymetry=(sur_siz_u-sur_siz_l)./(sur_siz_u+sur_siz_l);
            
            center_size=sum(center_mask,'all');
            center_surround_com_orientation=RFprops{n}.orientation_centers_of_mass;
            center_surround_distance=RFprops{n}.distance_btw_centers_of_mass;
        else
            rel_surround_strength=nan;
            center_size=nan;
            center_surround_com_orientation=nan;
            center_surround_distance=nan;
            vertical_surround_assymetry=nan;
            vertical_surround_size_assymetry=nan;
        end
        RF_param(n,18:23)=cat(2,center_size,rel_surround_strength,vertical_surround_assymetry,center_surround_distance,center_surround_com_orientation,vertical_surround_size_assymetry);
        
        waitbar(n/size(RF_param,1),w)
    end
    close(w)
    %
    RFa_bgA=cat(1,RFa_bgA,RFa_bg);
    RFProps_a=cat(1,RFProps_a,RFprops);
    RFa_param{rf}=RF_param;
    
    
end
RFa_bg=RFa_bgA;
sac_deg_tuning=sac_deg_tuningA;
clear RFa_bgA
RFa_param_key=RF_param_key;
RFa_param=cat(1,RFa_param{:});
%% make retinotopic fits
SNR = 10*log10(RFa_param(:,5)./RFa_param(:,6));
RFa_param(:,26:28)=nan;
RFa_param_key{26}='retinotopic_predEl';
RFa_param_key{27}='retinotopic_predAz';
RFa_param_key{28}='retinotopic_distance';
RFa_param_key{29}='SNR';
RFa_param(:,29)=SNR;
foptions = fitoptions('poly22');

for rf=1:numel(resultfiles)
    
    temp=SNR(RFa_param(:,9)==rf);
    if numel(temp)>6
        std_thresh=prctile(temp,85);
        
        quadr_intensity=1e-4; %curvature of plane
        foptions.Weights=SNR(RFa_param(:,9)==rf)-std_thresh;%use SNR as weight
        foptions.Weights(foptions.Weights<0)=0;%remove low SNR from fit
        foptions.Exclude=SNR(RFa_param(:,9)==rf)<std_thresh;
        foptions.Lower=[-inf -inf -inf -quadr_intensity -quadr_intensity -quadr_intensity];
        foptions.Upper=[inf inf inf quadr_intensity quadr_intensity quadr_intensity];
        
        [sf_az,gof1] = fit(double(RFa_param(RFa_param(:,9)==rf,11:12)),double(RFa_param(RFa_param(:,9)==rf,4)),'poly22',foptions);
        [sf_el,gof2] = fit(double(RFa_param(RFa_param(:,9)==rf,11:12)),double(RFa_param(RFa_param(:,9)==rf,3)),'poly22',foptions);
        
        ElAz_pred=cat(2,sf_el(RFa_param(RFa_param(:,9)==rf,11:12)),sf_az(RFa_param(RFa_param(:,9)==rf,11:12)));%predict retinotopic position of all ROIs
        
        retinotopic_distance=RFa_param(RFa_param(:,9)==rf,3:4)-ElAz_pred; 
        rf_bar_offd=sqrt(sum(retinotopic_distance.^2,2));
        RFa_param(RFa_param(:,9)==rf,26:28)=[ElAz_pred rf_bar_offd];
    end
end

RFa_param(:,30:31)=RFa_param(:,3:4);
RFa_param_key(30:31)={'pos_deg_el_notworldrotated','pos_deg_az_notworldrotated'};
%% make 1D RFs
peakTime2D=paramsA.relativeTime(RFa_param(:,8));
RF1d=cellfun(@(x) permute(mean(x(:,21:29,:,:),2),[1 4 2 3]),RFa_bg,'UniformOutput',0);
RF1d=cat(3,RF1d{:});
RF1d=RF1d-median(RF1d(:,paramsA.relativeTime>0,:),[1 2]);
RF1d=RF1d./max(abs(RF1d),[],[1 2]);
[c_peak_val, c_peaktemp] = max(abs(RF1d(:,paramsA.relativeTime<0,:)),[],[1 2],'linear');
[c_peakEl,c_peaktime,peakN]=ind2sub(size(RF1d(:,paramsA.relativeTime<inf,:)),c_peaktemp(:));
temp2=permute(var(RF1d(:,paramsA.relativeTime<inf,:),0,1),[3 2 1]);
[~, c_peaktimev] = max(temp2,[],2);
temp2=permute(var(RF1d(:,paramsA.relativeTime<inf,:),0,2),[3 1 2]);
[~, c_peakEl] = max(temp2,[],2);
peakTime1D=paramsA.relativeTime(c_peaktimev);
SNR1d=nan(size(c_peakEl));
peak_val = c_peak_val.^2;
for n=1:size(RF1d,3)
    noisevals = RF1d(:,c_peaktimev(n),n);
    nsel=c_peakEl(n)-15:c_peakEl(n)+15;
    nsel(nsel<1)=[];nsel(nsel>numel(noisevals))=[];
    pad_area=find(all(RF1d(:,:,n)==0,2));
    nsel=unique(cat(2,nsel,pad_area'));
    Ns=numel(noisevals)-numel(nsel);
    noisevals(nsel)=[];
    noise_var=sum(noisevals.^2)/Ns;
    SNR1d(n) = 10*log10(peak_val(n)/noise_var);
end

%% fit parameterizations to 1D RFs
RF_fittype_1d =  fittype('amp1*exp(-((x-x01).^2/(2*wx1^2)))-amp2*exp(-((x-x02).^2/(2*wx2^2)))', 'dependent',{'z'},'independent',{'x'},...
    'coefficients',{'amp1','x01','wx1','amp2','x02','wx2'});
fitparams_1d=nan(size(RF1d,3),7);
zero_t=find(paramsA.relativeTime>0,1);
steps=unique(min(0:24:size(RF1d,3)+400,size(RF1d,3)))';
w=waitbar(0,'fitting 1D RFs');
elvals=(1:size(RF1d,1))';
elvals2=linspace(1,size(RF1d,1),10*size(RF1d,1))';
for st=1:numel(steps)-1
    startst=steps(st)+1;
    endst=steps(st+1);
    parfor n=startst:endst
        if SNR1d(n)>10 && peakTime1D(n)<=0.2 && peakTime1D(n)>=-.7
            cFRF=double(RF1d(:,:,n));
            test=cFRF(:,c_peaktimev(n));%select peak time for fit
            [test_max,test_maxy]=max(test);
            [test_min,test_miny]=max(-test);
            off_cell=test_min-test_max>0;
            if off_cell %invert
                test2=-test;
                cFRF=-cFRF;
            else
                test2=test;
            end
            
            cFRF=cFRF-median(test2);
            test2=test2-median(test2);

            % start points for fit
            [pH,pL,pW]=findpeaks(test2,'NPeaks',1,'MinPeakHeight',max(test2)*.9);
            if ~isempty(pH)
                pW=2*sqrt(pW);
                StartPoint=[pH,     pL,     pW,     0,   pL,    10*pW];
                Lower=     [pH/2,     max(5,pL-20),  .25*pW, 0,   0,     pW];
                Upper=     [inf,   min(pL+20,numel(test2)-5),  inf,   2*pH,  numel(test2),   inf];
            else
                [pH,pL]=max(test2);
                StartPoint=[pH,     pL,     5,     0,   pL,    10];
                Lower=     [pH/2,     5,  1, 0,   0,     1];
                Upper=     [inf,   numel(test2)-5,  inf,   pH,  numel(test2),   inf];
                
            end
            RF_fitoptions=fitoptions(RF_fittype_1d);
            RF_fitoptions.StartPoint=StartPoint;
            RF_fitoptions.Lower= Lower;
            RF_fitoptions.Upper= Upper;
            [RF_fitresult, RF_gof] = fit( (1:size(test,1))', test2, RF_fittype_1d, RF_fitoptions);
            
            %get 1D parameterizations
            model=RF_fitresult(elvals);
            centvals=test2.*(model>0);
            survals=test2.*(model<0);
            com=[sum(elvals.*centvals)/sum(centvals) sum(elvals.*survals)/sum(survals)];
            tempamps=nan(size(cFRF,2),2);
            
            half_size=floor(numel(elvals)/2);
            center_mask=false(size(test));
            center_mask(elvals>=RF_fitresult.x01-2*RF_fitresult.wx1 & elvals<=RF_fitresult.x01+2*RF_fitresult.wx1)=true;
            center_mask2=false(size(elvals2));
            center_mask2(elvals2>=RF_fitresult.x01-2*RF_fitresult.wx1 & elvals2<=RF_fitresult.x01+2*RF_fitresult.wx1)=true;
            
            surround_mask=false(size(test));
            surround_mask(elvals>=RF_fitresult.x02-2*RF_fitresult.wx2 & elvals<=RF_fitresult.x02+2*RF_fitresult.wx2)=true;
            surround_mask(center_mask)=false;
            surround_mask_upper=surround_mask;
            
            surround_mask_upper(1:half_size)=false;
            surround_mask_lower=surround_mask;
            surround_mask_lower(end-half_size:end)=false;
            
            surround_mask_upper2=surround_mask;
            surround_mask_upper2(1:floor(RF_fitresult.x01))=false;
            surround_mask_lower2=surround_mask;
            surround_mask_lower2(ceil(RF_fitresult.x01):end)=false;
            sur_sum=abs(sum(test2(surround_mask)));
            cent_sum=abs(sum(test2(center_mask)));
            rel_surround_strength=abs(sur_sum./cent_sum);
            sur_sum_l=abs(sum(test2(surround_mask_lower)));
            sur_sum_u=abs(sum(test2(surround_mask_upper)));
            sur_sum_l2=abs(sum(test2(surround_mask_lower2)));
            sur_sum_u2=abs(sum(test2(surround_mask_upper2)));
            vertical_surround_assymetry=(sur_sum_u-sur_sum_l)./(sur_sum_u+sur_sum_l);
            vertical_surround_assymetry2=(sur_sum_u2-sur_sum_l2)./(sur_sum_u2+sur_sum_l2);
            center_size2=sum(center_mask2)/10;
            center_size=2*RF_fitresult.wx1;%^2
            center_surround_distance=diff(com);
            sur_siz_l=sum(surround_mask_lower,'all');
            sur_siz_u=sum(surround_mask_upper,'all');
            vertical_surround_size_assymetry=(sur_siz_u-sur_siz_l)./(sur_siz_u+sur_siz_l);
            
            %% collect 1d parameterization data
            if off_cell;off_cell_fac=-1;else;off_cell_fac=1;end
            Fits1d(n).fitresult=RF_fitresult;
            Fits1d(n).fitgof=RF_gof;
            Fits1d(n).ampst=off_cell_fac.*[1 -1].*tempamps; %flip on cells
            Fits1d(n).coefs=[off_cell_fac*RF_fitresult.amp1, RF_fitresult.x01,RF_fitresult.wx1,-off_cell_fac*RF_fitresult.amp2,RF_fitresult.x02,RF_fitresult.wx2];
            Fits1d(n).com=com;
            
            fitparams_1d(n,:)=cat(2,center_size,rel_surround_strength,vertical_surround_assymetry,center_surround_distance,sur_sum,center_size2,vertical_surround_assymetry2);
            
        else
            Fits1d(n).ampst=nan(size(RF1d(:,:,n),2),2);
            Fits1d(n).coefs=nan(1,6);
            Fits1d(n).com=nan(1,2);
        end
    end
    waitbar(st/(numel(steps)-1),w)
end
close(w)

clearvars -except c_peaktimev RF1d peakTime1D peakTime2D SNR1d fitparams_1d Fits1d RFProps_a RFa_param ch_siz_rate RFa_param_key RFa_bg resultfiles GeneralParams paramsA

%% save intermediate full collection (optional)
save(collected_savename,'c_peaktimev','RF1d','peakTime1D','peakTime2D','SNR1d','fitparams_1d','Fits1d','RFProps_a','RFa_param','ch_siz_rate','RFa_param_key','RFa_bg','resultfiles','GeneralParams','paramsA','-v7.3');
%% load full collection (optional)
load(collected_savename);

%%
%% inclusion criteria

axis_offset=0; %azimuth axis offset for world rotation
rotation=-20; % world rotation to compensate
min_R2=0.0; % minimal R2 fit with RF (unused)
min_SNR=15; % minimal SNR of RF
max_retinotopic_dist=20; %maximal distance (deg) of RF center to retinotopic position
t_peak_win=[-.6 0.1]; % accepted window of temporal peak of RF
max_rss=inf; %maximum relative surround strength (unused)
elmax=35; %maximal elevation (deg) of RF center to avoid edge
elmin=-35;  %minimal elevation (deg) of RF center to avoid edge


peakTime2D=peakTime2D(:);
peakTime1D=peakTime1D(:);
RF1d=cellfun(@(x) permute(mean(x(:,18:32,:,:),2),[1 4 2 3]),RFa_bg,'UniformOutput',0);
RF1d=cat(3,RF1d{:});
RF1d=RF1d-median(RF1d(:,paramsA.relativeTime>0,:),[1 2]);

SNR_1d=10*log10(squeeze(max(abs(RF1d),[],[1 2])./var(RF1d(:,paramsA.relativeTime>0.2 | paramsA.relativeTime<-0.7,:),0,[1 2])));
RF1d=RF1d./max(abs(RF1d),[],[1 2]);

SNR_2d=10*log10(RFa_param(:,5)./RFa_param(:,6));
%
if min_R2>0
fit_r2=nan(numel(Fits1d),2);
for n=1:numel(Fits1d)
    if ~isempty( Fits1d(n).fitgof)
        fit_r2(n,1)=Fits1d(n).fitgof.rsquare;
        fit_r2(n,2)=Fits1d(n).fitgof.rmse;
    end
    if RFProps_a{n}.center_exist==1
        R2vals=mk_R2_RF(RFProps_a{n});
        fit_r2(n,3)=max(R2vals.R2_cs,R2vals.R2_c);
        fit_r2(n,4)=R2vals.rmse_cs;
    end
    
end
end
%% rotate world

fprintf('rotation: %.0f\n',rotation);
R = makehgtform('yrotate',deg2rad(-rotation));
R = R(1:3,1:3);
coords=nan(size(RFa_param,1),3);
[coords(:,1),coords(:,2),coords(:,3)] = sph2cart(deg2rad(RFa_param(:,31)+axis_offset),deg2rad(RFa_param(:,30)),40) ;
rotcoords_cart=coords;
for c=1:size(coords,1)
    rotcoords_cart(c,:) = (R*coords(c,:)')';
end
rotcoords_sph=rotcoords_cart;
[rotcoords_sph(:,1),rotcoords_sph(:,2),rotcoords_sph(:,3)]=cart2sph(rotcoords_cart(:,1),rotcoords_cart(:,2),rotcoords_cart(:,3));
RFa_param(:,3:4)=rad2deg(rotcoords_sph(:,[2 1]))-[0 axis_offset];

%% get selection




selNpos_off=RFa_param(:,28)<max_retinotopic_dist;
selN_SNR_2d=SNR_2d>min_SNR;
selN_SNR_1d=SNR1d>min_SNR;
selN_R2_2d=fit_r2(:,3)>=min_R2;
selN_R2_1d=fit_r2(:,1)>=min_R2;
sel_2d_rss=RFa_param(:,19)<max_rss;
sel_1d_rss=fitparams_1d(:,2)<max_rss;
sel_fit_2d=~any(isnan(RFa_param(:,18:20)),2);
sel_fit_1d=~any(isnan(fitparams_1d),2);
sel_el=RFa_param(:,30)<=elmax & RFa_param(:,30)>=elmin;
peakpos=cat(1,Fits1d(:).coefs);
peakpos=peakpos(:,[2 5]);
sel_sur_peak_off=peakpos(:,2)>0 & peakpos(:,2)<=51; %surround peak within cropped RF
selNtpeak_2d=(peakTime2D>=t_peak_win(1) & peakTime2D<=t_peak_win(2));
selNtpeak_1d=(peakTime1D>=t_peak_win(1) & peakTime1D<=t_peak_win(2));

%selections
subsel_2d=selN_SNR_2d & selNpos_off & selNtpeak_2d & sel_el;
subsel_1d=selN_SNR_1d & selNpos_off & selNtpeak_1d & sel_el;
subsel_2d_param=sel_fit_2d & selN_SNR_2d & selNpos_off & selNtpeak_2d & selN_R2_2d & sel_2d_rss & sel_el;
subsel_1d_param=sel_sur_peak_off & sel_fit_1d & selN_SNR_1d & selNpos_off & selNtpeak_1d & selN_R2_1d & sel_1d_rss & sel_el;
csubsel=subsel_1d_param; %used selection
fprintf('selected %.3f (2D), %.3f (2D param), %.3f (1D), %.3f (1D param)\n',mean(subsel_2d),mean(subsel_2d_param),mean(subsel_1d),mean(subsel_1d_param))



%% select RFs & parameters, invert positive center RFs
RF=RFa_bg(csubsel);
RF1=RF1d(:,:,csubsel);
peakTp=c_peaktimev(csubsel);

RF_parameters_2D=RFa_param(csubsel,:);
RF_parameters_1D=fitparams_1d(csubsel,:);
RF_2D=cell(size(RF));
RF_1D=cell(size(RF));
RF_on_cell_2D=false(size(RF));
RF_on_cell_1D=false(size(RF));
for n=1:numel(RF)
    RF_2D{n}=RF{n}(:,:,:,peakTp(n));%2D RF at peak time
    
    RF_1D{n}=RF1(:,peakTp(n),n);%1D RF at peak time
    %invert positive center RFs
    if mean(RF_2D{n}(26-1:26+1,26-1:26+1),'all')>0 %center_mag(n)>0 %
        RF_2D{n}=-RF_2D{n};
        RF_on_cell_2D(n)=true;
    end
    if mean(RF_1D{n}(26-1:26+1))>0 %center_mag(n)>0
        RF_1D{n}=-RF_1D{n};
        RF_on_cell_1D(n)=true;
    end
end
RF_2D=cat(3,RF_2D{:});
RF_1D=cat(3,RF_1D{:});
RF_SNR_1D=SNR1d(csubsel);
%% save selected file
RF_parameters_key=RFa_param_key(1:31);
save(selected_savename,'saccade_axis','saccade_parameters','saccade_params_key','RF_SNR_1D','RF_parameters_2D','RF_parameters_1D','RF_2D','RF_1D','RF_on_cell_1D','RF_on_cell_2D','RF_parameters_key' )