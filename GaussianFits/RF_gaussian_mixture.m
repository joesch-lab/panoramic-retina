%% use this script to fit a mixture of two gaussions to RF.
% RF of one neuron is the 3D array [rows,columns, time_samples]
% two gaussians (with positive and negative amplitudes) will be fitted
% then at every time step the mixture
% coefficients for these gaussins will be fitted

% General description:
% First I clean the RF by setting to 0 (and 1) those pixels that have weak 
% (strong) correlation with the peak in the RF (threshold param corr_threshold).
% Then I fit gaussians to the median of the RF around the peak (time_peak +/- time_window).
%  (1) fit one 2D gaussian G1 to the median_RF
%  (2) fit the sum of gaussians, G1+G2
%     a) only the amplitude of G1 can be adjusted, other parameters stay fixed.
%     b) the amplitudes G2 and G1 are of opposite signs sign(A(G2))=(-1)*sign(A(G1))
%     c) there are two cases for fitting the 2nd gaussian: 
%          - G2 can be the surround of G1 (in this case abs(A(G2))<<abs(A(G1))), 
%            in this case the center of G2 is restricted to be within 2*r from the center of G1;
%          - G2 can be an independent component (when A(G2) ==A(G1)
%            and centers of G1 and G2 are well separated).
%            in this case there are no restrictions on the location of the center of G2.
%  (3) for each time step, fit scaling parameters of the gaussian mixture lambda1 and lambda2; 
%      I allow both lambdas to be in the interval [-Inf,Inf], however sign(lamda1) = -sign(lamda2). 
%      In this way I can model the ON-OFF or OFF-ON dynamics of the center 
%      while enforcing the second component to be always of the opposite sign.
%  (4) using the widths of the gaussians I create a mask of the center and of the surround. 
%      I compute and store the correlations between the mean values of the data 
%      convolved by masks across the time. This correlation value might indicate 
%      the goodness of fit and could be used to filter out false-positives.

% O.Symonova
% 11.02.2021

function RFprops = RF_gaussian_mixture(rf, resfolder)

     
    %if RF array has different sequence shift dimentions accordingly
    cluster_info_tsv_file=[];
%     fileRF='C:\DATA\Data_for_Olga\retina\rfs_of_20201027-134944_using_spks.mat';
%     [datafolder,filename,fileext]=fileparts(fileRF);
    
    if nargin > 1
        %folder where to store results
    %     resfolder=fullfile(datafolder,'test');
%         resfolder = "~/codes/retina-analysis/gaussian-mixture/res";
        %make a figure of the model fit for each cluster
        make_figure=1;
    else
        make_figure = 0;
    end
   
    
%     load(fileRF);
if isfield(rf, 'bg_mean')
    RF=rf.filter-rf.bg_mean;
else
    RF = rf.filter;
end
    RF=squeeze(RF);
    addpath("../RFAnalysis/");
    RF = normalize_filter(RF);

    %expecteds sequence of dimensions is [nrow, ncol, ntframes, nclusters]
    RF=shiftdim(RF,1);
    RF=double(RF);

    
    %the method requires the mean of the RF to be at zero (for curve fitting to work)
    %if RF has range [0,1] or [0,255] then set substract_mean=1 and the
    %method will transform it to [-0.5, 0.5] or [-127,127] by substracting
    %the the mean.
    substract_mean=1;
    
    [nrow,ncol,nhist,ncl]=size(RF);
    nrowncol=nrow*ncol;

    %params for for the gaussian fit
   
    %init guess of the size of RF is 10 degrees, std is half width    
    rfw0=15;%round((608*10/220)/2);  
    % limits of the max and minimum width of a RF
    %minRF is 1 degree
    minwRF=3;%(608*1/220)/2;
    %maxRF is 110 degree
    maxwRF=25;%round((608*110/220)/2);
    %how close to the border RF is allowed
    border_proximity = 0;%100;
    %time window around the peak to compute the image, where to fit
    %gaussians
    twin=3;
    %time window (tpeak-twin_corr;tpeak+twin_corr) around the peak to crate the correlation mask
    % if the window is too small a lot of pixels might appear uncorrelated;
    % if it is too large random fluctuations before and after RF might
    % dominate; 
    twin_corr=10;
    
    %correlation threshold, pixels with correlation below this value will
    %not be used in the fitting
    corr_threshold = 0.4;
    
    %RF variance should be a factor above the mean variance, the higher
    %this threshold, the more reliable are the results (fewer false positives but more false negatives)   
    min_aboveavvar=5;
    
    %the relative strenth of positve and negative components; if amplitude
    %of one component is much smaller then the other, then it is
    %center_surround configuration
    min_rel_comp_strength = 0.75;
    
    %set verbous mode for curve fitting off
    myopts = optimoptions(@lsqcurvefit,'display','off','UseParallel',false);
           
    %grid space for the whole image
    [Xt Yt] = meshgrid(1:ncol,1:nrow);
    XYt=[];
    XYt(:,1) = Xt(:)'; %convert image into column format
    XYt(:,2) = Yt(:)';
    %square format of the image grid
    XYtim=[];
    XYtim(:,:,1)=Xt;
    XYtim(:,:,2)=Yt;
    
    %init frame pixels to better estimate noise
    frame=ones(nrow,ncol);
    frame_width = 5;
    frame(frame_width+1:nrow-frame_width+1,frame_width+1:ncol-frame_width+1)=0;
    frame_size=nnz(frame);
    peak2noise_log_threshold=15;
     
    %init properties structure
% %     RFprops=[];
%    selected = [98 97 95 94 91 87 85 83 82 81 80 79 74 73 72 70 67 66 62 61 59 58 56 54 52 51 50 47 45 44 43 41 40 37 36 35 32 31 29 27 ];
%    selected = [96 107 131 177 187];
    for i=1:ncl
%        for sel = 1:numel(selected)
%            i = selected(sel);
%     for i=2168:2169
        disp(['Processing cluster ',num2str(i),'.']);
        RFi=RF(:,:,:,i);
        
        RFprops(i).error='';
        
        
        %find the variance of the RF
        %if max var is low compared to the mean var, skip the cluster
        vari=var(RFi,0,3);
        varinnz=numel(vari(vari==0));
        if varinnz/nrowncol>0.1
            RFprops(i).error='Location of the patch is at the border.';
            continue;
        end
        meanv=mean(vari(:));
        [maxv,maxind]=max(vari(:));
        

        
        
        RFprops(i).mean_variance = meanv;
        RFprops(i).max_variance = maxv;
        RFprops(i).max_mean_variance = maxv/meanv;
%         if maxv<min_aboveavvar*meanv
%             RFprops(i).error='low variance of RF. No gaussians initialized.';
%             continue;
%         end
        if isnan(maxv)
            RFprops(i).error='RF contains nan values.';
            continue;
        end
        
        %else init centers of both gaussians at the peak of the variance
        [ypeak, xpeak]=ind2sub([nrow,ncol],maxind);
        %if the peak is close to the border skip the component
        if xpeak <  border_proximity || xpeak > ncol-border_proximity || ypeak < border_proximity || ypeak > nrow-border_proximity
            RFprops(i).error='peak close to border';
            continue;
        end
            
        
        %find the time when there is most variance in the background
        vart=zeros(1,nhist);
        for hi=1:nhist
            rfi_t = RFi(:,:,hi);
            vart(hi)=var(rfi_t(:));
        end
        [maxvart, tpeak]=max(vart);
        
        if substract_mean
            RFi_mean=median(RFi,3); 
            RFi_nomean = RFi-RFi_mean;
        else
            RFi_nomean = RFi;
        end
        
        %compute new Peak-to-Noise ratio
        %noise is estimated as the variance of the pixels in a thin outer border
        %temporal trace of the peak
        cij=RFi_nomean(ypeak,xpeak,:);

        %when does the peak happen
        [c_peak_val, c_peaktime] = max(abs(cij));        
        peak_val = c_peak_val^2;        
        noisevals = RFi_nomean(:,:,c_peaktime);
        noisevals = noisevals(frame>0);
        noise_var=sum(noisevals.^2)/frame_size;
        RFprops(i).peak2noise = peak_val/noise_var;
        RFprops(i).peak2noise_log = 10*log10(RFprops(i).peak2noise);
        if RFprops(i).peak2noise_log<peak2noise_log_threshold %too noisy
           RFprops(i).error='too noisy';
            continue;
        end
        
        
        %compute correlations with the peak pixel
        tst=max(1,tpeak-twin_corr);
        ten=min(nhist,tpeak+twin_corr);        
        cij=RFi_nomean(ypeak,xpeak,tst:ten);
        corrim=zeros(nrow,ncol);
        for ii=1:nrow
            for jj=1:ncol                
                vij=RFi_nomean(ii,jj,tst:ten);
                if any(vij)                    
                    cr=corrcoef(vij(:),cij(:));
                    corrim(ii,jj)=cr(1,2);  
                end
            end
        end
        %correlation mask
        corrmask=corrim;
        corrmask(abs(corrmask)<corr_threshold)=0;
        corrmask(corrmask~=0)=1;
        RFprops(i).corr_mask = corrmask;
        corrmask_colfmt=logical(corrmask(:)); %correlation mask in column format
        
        %clean RF array, use pixels which are correlated with the peak
        RFi_nomean_clean=RFi_nomean.*corrmask;
        
        %find max and min amplitude at the peak and its location
        [RFijmax,tmaxind] = max(RFi_nomean_clean(:));
        [RFijmin,tminind] = min(RFi_nomean_clean(:));
        [ypeakp, xpeakp, tpeakp] = ind2sub(size(RFi_nomean_clean),tmaxind);            
        [ypeakn, xpeakn, tpeakn] = ind2sub(size(RFi_nomean_clean),tminind);            
        
        %first fit one gaussian around the biggest peak
        if RFijmax>abs(RFijmin) %ON cell
            ypeak=ypeakp; xpeak=xpeakp; tpeak=tpeakp;
            Acenter=RFijmax;            
            Bminpos_c=[0,border_proximity,border_proximity,minwRF,minwRF,0];
            Bmaxpos_c =[Inf,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];
        else %OFF cell
            ypeak=ypeakn; xpeak=xpeakn; tpeak=tpeakn;
            Acenter=RFijmin;
            Bminpos_c = [-Inf,border_proximity,border_proximity,minwRF,minwRF,0];
            Bmaxpos_c = [0,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];        
        end
        
        tst=max(1,tpeak-twin);
        ten=min(nhist,tpeak+twin);        
        RFaroundpeak=median(RFi_nomean(:,:,tst:ten),3); 
        RFprops(i).RFaroundpeak = RFaroundpeak;
%         RFaroundpeak_clean=RFaroundpeak.*corrmask;
        RFaroundpeak_clean=RFaroundpeak(:); 
        RFaroundpeak_clean=RFaroundpeak_clean(corrmask_colfmt); %data to fit, only nonzero values
        XYt_i=XYt(corrmask_colfmt,:); %coordinates of nonzero values
        
%         data = RFaroundpeak;
%         maxamp=max(max(data(:)),abs(min(data(:))));
%         datai=data./maxamp;      % in the range[-1,1] 
%         imagesc(datai, [-1,1]);
%         colormap(bluered);  
%         disp(['Max_var/mean_var:',num2str(maxv/meanv),'.']);
%         input("");       
%         continue;
%         
        
        %fit center
        B0c =  [Acenter, xpeak,ypeak, rfw0,rfw0, 0.0];
        %[xc,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian_Rot_center,B0c,XYt,RFaroundpeak_clean,Bminpos_c,Bmaxpos_c, myopts); 
        [xc,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian_Rot_center_colfmt,B0c,XYt_i,RFaroundpeak_clean,Bminpos_c,Bmaxpos_c, myopts); 
                  

        %check if it is center-surround or two independent conponents 
        case_two_strong_components =min(RFijmax,abs(RFijmin))/max(RFijmax,abs(RFijmin)) > min_rel_comp_strength;
        d=sqrt((ypeakp-ypeakn)^2+(xpeakp-xpeakn)^2);
        rc=2*min(xc(4:5));
        %distance btw pos and neg peaks is further than 2std of the center
        %gaussian
        case_separated_peaks = d>rc;
        
        if case_two_strong_components && case_separated_peaks
            disp(['Cluster ',num2str(i),': two independent components (',num2str(min(RFijmax,abs(RFijmin))/max(RFijmax,abs(RFijmin))),').']);
            %fit two independent gaussins (the center of the smaller is not locked to the center of the bigger)
            %fit the 2nd gaussiand and allow to re-fit the amplitude of the first
            if xc(1)>0 %on-cell                           
                Bminpos=[0,xc(2:6),-Inf,border_proximity,border_proximity,minwRF,minwRF,0];
                Bmaxpos =[Inf,xc(2:6),0,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];                
                Asurround = RFijmin;
                xs=xpeakn; ys = ypeakn;
            else %off-cell
                Bminpos=[-Inf,xc(2:6),0,border_proximity,border_proximity,minwRF,minwRF,0];
                Bmaxpos =[0,xc(2:6),Inf,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];
                Asurround = RFijmax;
                xs=xpeakp; ys = ypeakp;
            end
            B0=[xc(1:6),Asurround, xs,ys, rfw0,rfw0, 0.0];
        else %fit two gaussins (the center of the smaller should be with the 2std of the bigger)
            disp(['Cluster ',num2str(i),': center-surround']);
            xs_min = max(1+border_proximity,xc(2)-rc);
            xs_max = min(ncol-border_proximity, xc(2)+rc);
            ys_min = max(1+border_proximity,xc(3)-rc); 
            ys_max = min(nrow-border_proximity,xc(3)+rc);
            if xc(1)>0 %on-cell
                Bminpos=[0,xc(2:6),-Inf,xs_min,ys_min,xc(4:5),0];
                Bmaxpos =[Inf,xc(2:6),0,xs_max,ys_max,maxwRF,maxwRF,2*pi];                
                Asurround = RFijmin;                
            else %off-cell
                Bminpos=[-Inf,xc(2:6),0,xs_min,ys_min,xc(4:5),0];
                Bmaxpos =[0,xc(2:6),Inf,xs_max,ys_max,maxwRF,maxwRF,2*pi];
                Asurround = RFijmax;                
            end
            B0=[xc(1:6),Asurround, xc(2:3),xc(4:5)*2, 0.0];            
        end
        
        %uncomment for the faster run, might not find the global optimum
        % [x,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian_Sum_Rot_center_surround_sep,B0,XYt,RFaroundpeak_clean,Bminpos,Bmaxpos, myopts);

        %use the option below to seed the optimisation problem at several
        %locations, might take awhile
        problem = createOptimProblem('lsqcurvefit',...
                'x0',B0,'objective',@Gaussian_Sum_Rot_center_surround_sep_colfmt,'xdata',XYt_i,'ydata',RFaroundpeak_clean,'lb',Bminpos,'ub',Bmaxpos,'options',myopts);
        ms = MultiStart;
        [x,f] = run(ms,problem,20);       

        RFprops(i).center=x(1:6);
        RFprops(i).surround=x(7:end);
        center=RFprops(i).center;
        surround=RFprops(i).surround;
        
        %compute correlation of the pixels in the surround with the pixels
        %in the center; %this value might indicate the goodness of fit
        center_mask=make_ellipse_mask(nrow,ncol, center(2:end));
        surround_mask=make_ellipse_mask(nrow,ncol, surround(2:end));
        surround_mask=surround_mask-center_mask;
        surround_mask(surround_mask<0)=0;
        RFprops(i).center_mask = center_mask;
        RFprops(i).surround_mask = surround_mask;        

        mean_center=zeros(1,nhist);
        mean_surround=zeros(1,nhist);
        for ii=1:nhist
            ci=RFi_nomean(:,:,ii).*center_mask;
            si=RFi_nomean(:,:,ii).*surround_mask;       
            mean_center(ii)=mean(nonzeros(ci));
            mean_surround(ii)=mean(nonzeros(si));            
        end
        tst=max(1,tpeak-twin_corr);
        ten=min(nhist,tpeak+twin_corr);
        cr=corrcoef(mean_center(tst:ten),mean_surround(tst:ten));
        RFprops(i).surround_center_corr = cr(1,2);
        
        model_center_surround = Gaussian_Sum_Rot_center_surround_sep([center, surround],XYtim);
        RFprops(i).model = model_center_surround;

        if make_figure
            filename=['n',num2str(i),'_tpeak',num2str(tpeak)]; 
            figure_name = fullfile(resfolder,filename);
            title_str=['Cluster ',num2str(i),', c-s corr ',num2str(cr(1,2)),'.'];
            make_data_model_figure(RFaroundpeak,model_center_surround, center, surround,title_str,figure_name);
        end
        
        %every time step fit scale parameters of the gaussians
        lambdas=[0.5,0.5];
        lambdas_min=[-Inf,-Inf];
        lambdas_max=[Inf,Inf];
        center_dynamics=zeros(nhist,1);
        surr_dynamics=zeros(nhist,1); 
        datamix.param = [RFprops(i).center, RFprops(i).surround];
        datamix.xdata = XYt_i;
        for ti=1:nhist
            RFii = RFi_nomean(:,:,ti);            
            RFii = RFii(:);
            RFii =RFii(corrmask_colfmt); %data to fit, only nonzero values
            [lambdas_fit,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian_Mixture_Rot_oppositesign_colfmt,lambdas,datamix,RFii,lambdas_min,lambdas_max, myopts); 
            center_dynamics(ti)=lambdas_fit(1);
            surr_dynamics(ti)=(-1)*sign(lambdas_fit(1))*lambdas_fit(2);        
        end
        
        RFprops(i).center_dynamics=center_dynamics;
        RFprops(i).surround_dynamics=surr_dynamics;
    end
    
if ~isempty(cluster_info_tsv_file)
    RFprops = set_depth_property(RFprops, cluster_info_tsv_file);    
end
% propsfile=fullfile(resfolder,'RF_props.mat');
% save(propsfile, 'RFprops');

RFprops = gaussian_mixture_RF_properties(RF,RFprops);

end

function elmask= make_ellipse_mask(nrow,ncol, el)
    rp1=round(el(3)*2);
    rp2=round(el(4)*2);
    alpha_p=(el(5)*180)/pi;

    col_p=round(el(1));
    if col_p > ncol; col_p = ncol; end
    if col_p < 1; col_p = 1; end
    row_p=round(el(2));
    if row_p > nrow; row_p = nrow; end
    if row_p < 1; row_p = 1; end
    
    mystrel=strel('disk', max(rp1,rp2),8);
    immystrel=mystrel.Neighborhood;
    immystrel=imresize(immystrel,[2*rp2, 2*rp1],'nearest');
    immystrel=imrotate(immystrel,-alpha_p,'nearest');
    elmask=zeros(nrow,ncol);
    elmask(row_p,col_p)=1;
    elmask=imdilate(elmask,immystrel);
end

function make_data_model_figure(data, model, center, surround, title_str, filename)

    rb = bluered();
    fig= figure;
    suptitle(title_str);
    subplot(1,2,1);
        %find the largest amplitude
        maxamp=max(max(data(:)),abs(min(data(:))));
        datai=data./maxamp;      % in the range[-1,1] 
        imagesc(datai, [-1,1]);
        colormap(rb);        
        %plot ellipses
        hold on;         
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
        axis tight;
        axis image;

        

    subplot(1,2,2);
        %find the largest amplitude
        maxamp=max(max(model(:)),abs(min(model(:))));
        modeli=model./maxamp;      % in the range[-1,1] 
        imagesc(modeli, [-1,1]);
        axis image;
        colormap(rb);
        axis tight;        
        
    for ext = [".pdf"] % 
        outname = filename + ext;
        saveas(fig,outname);
    end
    close(fig);
end

function map = bluered()
m=256;
c=1;

if mod(m,2)
	z = [0,0,0];
	m2 = floor(m/2);
else
	z = zeros([0,3]);
	m2 = m/2;
end
map = [repmat([0,0,1],[m2,1]);z;repmat([1,0,0],[m2,1])];
r = repmat(abs(linspace(1,-1,m)).^c,[3,1])';
map = map.*r + 1 - r;
end


function RFprops = set_depth_property(RFprops, cluster_info_tsv_file)
    id_depth = get_clusters_depth(cluster_info_tsv_file);
    ncl=length(RFprops);
    ncl_depth = size(id_depth,1);
    if ncl~=ncl_depth
        warning(['The number of clusters in the info file and property',...
        'structure do not match. Could not assign depth property.']);
        return;
    end
    for i=1:ncl
        RFprops(i).id=id_depth(i,1);
        RFprops(i).depth=id_depth(i,2);
    end    
end




function [xn,yn] = rotate_pts(pt,ptcenter, alpha)
    pttemp=pt-ptcenter;
    xn = pttemp(1)*cos(alpha) - pttemp(2)*sin(alpha);
    yn = pttemp(1)*sin(alpha) + pttemp(2)*cos(alpha);
    xn=xn+ptcenter(1);
    yn=yn+ptcenter(2);
end
