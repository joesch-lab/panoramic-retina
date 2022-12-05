%% This script fits diff of Gaussian model to one spatial RF
%For detailed comments see RF_gaussian_mixture.m which fits gaussians to 3D
%RFs.

% O.Symonova, D. Gupta
% 11.02.2021

function RFprops = spatial_gprops(rf)

    
    %the method requires the mean of the RF to be at zero (for curve fitting to work)
    %if RF has range [0,1] or [0,255] then set substract_mean=1 and the
    %method will transform it to [-0.5, 0.5] or [-127,127] by substracting
    %the the mean.
    substract_mean=1;
    
%     [nrow,ncol,nhist,ncl]=size(RF);
    
    [nrow,ncol] = size(rf);
    
    
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
    
%     %time window around the peak to compute the image, where to fit
%     %gaussians
%     twin=3;
%     %time window (tpeak-twin_corr;tpeak+twin_corr) around the peak to crate the correlation mask
%     % if the window is too small a lot of pixels might appear uncorrelated;
%     % if it is too large random fluctuations before and after RF might
%     % dominate; 
%     twin_corr=10;
%     
%     %correlation threshold, pixels with correlation below this value will
%     %not be used in the fitting
%     corr_threshold = 0.4;
%     
%     %RF variance should be a factor above the mean variance, the higher
%     %this threshold, the more reliable are the results (fewer false positives but more false negatives)   
%     min_aboveavvar=5;
%     
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
     
    %init properties structure
% %     RFprops=[];
%    selected = [98 97 95 94 91 87 85 83 82 81 80 79 74 73 72 70 67 66 62 61 59 58 56 54 52 51 50 47 45 44 43 41 40 37 36 35 32 31 29 27 ];
%    selected = [96 107 131 177 187];
%        for sel = 1:numel(selected)
%            i = selected(sel);
%     for i=2168:2169
        
        RFprops.error='';
        
        
%         %find the variance of the RF
%         %if max var is low compared to the mean var, skip the cluster
%         vari=var(rf,0,3);
%         varinnz=numel(vari(vari==0));
%         if varinnz/nrowncol>0.1
%             RFprops.error='Location of the patch is at the border.';
%             continue;
%         end
%         meanv=mean(vari(:));
%         [maxv,maxind]=max(vari(:));
%         
% 
%         
%         
%         RFprops.mean_variance = meanv;
%         RFprops.max_variance = maxv;
%         RFprops.max_mean_variance = maxv/meanv;
%         if maxv<min_aboveavvar*meanv
%             RFprops.error='low variance of RF. No gaussians initialized.';
%             continue;
%         end
%         if isnan(maxv)
%             RFprops.error='RF contains nan values.';
%             continue;
%         end
%         
        %else init centers of both gaussians at the peak of the variance
%         [ypeak, xpeak]=ind2sub([nrow,ncol],maxind);
        
        ypeak = ceil(nrow/2);
        xpeak = ceil(ncol/2);
        
%         %if the peak is close to the border skip the component
%         if xpeak <  border_proximity || xpeak > ncol-border_proximity || ypeak < border_proximity || ypeak > nrow-border_proximity
%             RFprops.error='peak close to border';
%             continue;
%         end
        
        %find the time when there is most variance in the background
%         vart=zeros(1,nhist);
%         for hi=1:nhist
%             rfi_t = rf(:,:,hi);
%             vart(hi)=var(rfi_t(:));
%         end
%         [maxvart, tpeak]=max(vart);
        
%         if substract_mean
%             RFi_mean=median(rf,3); 
%             RFi_nomean = rf-RFi_mean;
%         else
%             RFi_nomean = rf;
%         end
        
%         %compute correlations with the peak pixel
%         tst=max(1,tpeak-twin_corr);
%         ten=min(nhist,tpeak+twin_corr);        
%         cij=RFi_nomean(ypeak,xpeak,tst:ten);
%         corrim=zeros(nrow,ncol);
%         for ii=1:nrow
%             for jj=1:ncol                
%                 vij=RFi_nomean(ii,jj,tst:ten);
%                 if any(vij)                    
%                     cr=corrcoef(vij(:),cij(:));
%                     corrim(ii,jj)=cr(1,2);  
%                 end
%             end
%         end
%         %correlation mask
%         corrmask=corrim;
%         corrmask(abs(corrmask)<corr_threshold)=0;
%         corrmask(corrmask~=0)=1;
%         RFprops.corr_mask = corrmask;
%         corrmask_colfmt=logical(corrmask(:)); %correlation mask in column format
        corrmask = true(size(rf));
        corrmask_colfmt = corrmask(:);
%         
        %clean RF array, use pixels which are correlated with the peak
%         RFi_nomean_clean=RFi_nomean.*corrmask;
        
        RFi_nomean_clean = rf;
        
        %find max and min amplitude at the peak and its location
        [RFijmax,tmaxind] = max(RFi_nomean_clean(:));
        [RFijmin,tminind] = min(RFi_nomean_clean(:));
        [ypeakp, xpeakp] = ind2sub(size(RFi_nomean_clean),tmaxind);            
        [ypeakn, xpeakn] = ind2sub(size(RFi_nomean_clean),tminind);            
        
        %first fit one gaussian around the biggest peak
        if RFijmax>abs(RFijmin) %ON cell
            ypeak=ypeakp; xpeak=xpeakp; 
            Acenter=RFijmax;            
            Bminpos_c=[0,border_proximity,border_proximity,minwRF,minwRF,0];
            Bmaxpos_c =[Inf,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];
        else %OFF cell
            ypeak=ypeakn; xpeak=xpeakn; 
            Acenter=RFijmin;
            Bminpos_c = [-Inf,border_proximity,border_proximity,minwRF,minwRF,0];
            Bmaxpos_c = [0,ncol-border_proximity,nrow-border_proximity,maxwRF,maxwRF,2*pi];        
        end
        
%         tst=max(1,tpeak-twin);
%         ten=min(nhist,tpeak+twin);        
%         RFaroundpeak=median(RFi_nomean(:,:,tst:ten),3); 
%         RFprops.RFaroundpeak = RFaroundpeak;
%         RFaroundpeak_clean=RFaroundpeak.*corrmask;
%         RFaroundpeak_clean=RFaroundpeak(:); 
%         RFaroundpeak_clean=RFaroundpeak_clean(corrmask_colfmt); %data to fit, only nonzero values
            RFaroundpeak_clean = rf(corrmask_colfmt);
        
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
            disp(['Cluster ',': center-surround']);
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

        RFprops.center=x(1:6);
        RFprops.surround=x(7:end);
        center=RFprops.center;
        surround=RFprops.surround;
        
        %compute correlation of the pixels in the surround with the pixels
        %in the center; %this value might indicate the goodness of fit
%         center_mask=make_ellipse_mask(nrow,ncol, center(2:end));
%         surround_mask=make_ellipse_mask(nrow,ncol, surround(2:end));
%         surround_mask=surround_mask-center_mask;
%         surround_mask(surround_mask<0)=0;
%         RFprops.center_mask = center_mask;
%         RFprops.surround_mask = surround_mask;        
% 
%         mean_center=zeros(1,nhist);
%         mean_surround=zeros(1,nhist);
%         for ii=1:nhist
%             ci=RFi_nomean(:,:,ii).*center_mask;
%             si=RFi_nomean(:,:,ii).*surround_mask;       
%             mean_center(ii)=mean(nonzeros(ci));
%             mean_surround(ii)=mean(nonzeros(si));            
%         end
%         tst=max(1,tpeak-twin_corr);
%         ten=min(nhist,tpeak+twin_corr);
%         cr=corrcoef(mean_center(tst:ten),mean_surround(tst:ten));
%         RFprops.surround_center_corr = cr(1,2);
        
        model_center_surround = Gaussian_Sum_Rot_center_surround_sep([center, surround],XYtim);
        RFprops.model = model_center_surround;

%         if make_figure
%             filename=['n',num2str(i),'_tpeak',num2str(tpeak)]; 
%             figure_name = fullfile(resfolder,filename);
%             title_str=['Cluster ',num2str(i),', c-s corr ',num2str(cr(1,2)),'.'];
%             make_data_model_figure(RFaroundpeak,model_center_surround, center, surround,title_str,figure_name);
%         end
        
        %every time step fit scale parameters of the gaussians
%         lambdas=[0.5,0.5];
%         lambdas_min=[-Inf,-Inf];
%         lambdas_max=[Inf,Inf];
%         center_dynamics=zeros(nhist,1);
%         surr_dynamics=zeros(nhist,1); 
%         datamix.param = [RFprops.center, RFprops.surround];
%         datamix.xdata = XYt_i;
%         for ti=1:nhist
%             RFii = RFi_nomean(:,:,ti);            
%             RFii = RFii(:);
%             RFii =RFii(corrmask_colfmt); %data to fit, only nonzero values
%             [lambdas_fit,resnorm,residual,exitflag] = lsqcurvefit(@Gaussian_Mixture_Rot_oppositesign_colfmt,lambdas,datamix,RFii,lambdas_min,lambdas_max, myopts); 
%             center_dynamics(ti)=lambdas_fit(1);
%             surr_dynamics(ti)=(-1)*sign(lambdas_fit(1))*lambdas_fit(2);        
%         end
%         
%         RFprops.center_dynamics=center_dynamics;
%         RFprops.surround_dynamics=surr_dynamics;
    
% propsfile=fullfile(resfolder,'RF_props.mat');
% save(propsfile, 'RFprops');

RFprops = extra_gprops(rf,RFprops);

end



function RFprops = extra_gprops(rf,RFprops)
%% this function takes the values of the guassians fitted to the RF
%  and computes different charachteristics using these gaussians. 
%  The new values are added as new fields to the RFprops structure.
%  The list of computed values:
%  - distance and angle between pos and neg guassuan peaks
%  - distance and angle between center of mass of pos and neg components
%  - time between the pos and neg peak
%  - half width of the center temporal dynamics
%  - masks of pos and neg gaussians (ellipses of 2sigma)
%  - average data and model values within the masks
%  - maximum data and model values within the masks
%  - average data and model values of only positive and only negative pixels 
%    within the masks
%  - number of positive and negative pixels within the masks (data and model)

% O.Symonova
% 11.02.2021

%%     %RFprops contains the following fields: 
%     fifilelename %name of the RF raw file
%     error %error message
% params of pos and negative guassians
%     center (the stronger component)
%     surround
% scaling params of the mixture
%     center_dynamics
%     center_dynamics
%     
%     correlation btw center and surround around the peak of the center     
%     
% 
%     
%     %parameters computed from fitting gaussian to the dynamics
%     %half width of temporal dynamics
%     hw_pos, hw_neg
%     %heigth and location of the peak of gaussians
%     peakheight_pos, peakheight_neg, peakloc_pos, peakloc_neg
%%
%     res_folder ='C:\DATA\Data_for_Olga\retina\';
    %the method requires the mean of the RF to be at zero (for curve fitting to work)
    %if RF has range [0,1] or [0,255] then set substract_mean=1 and the
    %method will transform it to [-0.5, 0.5] or [-127,127] by substracting
    %the the mean.
    substract_mean=1;    
    
    %set the flag to -1 for components with an error message
    RFprops = set_void_components(RFprops);
    
%     [nrow,ncol, nhist, ncl]=size(RF);
    [nrow, ncol] = size(rf);
    %create the grid to compute model values
    [X, Y] = meshgrid(1:ncol,1:nrow);
    Xt = X;
    Yt = Y;
    XY=[];
    XY(:,:,1) = X;
    XY(:,:,2) = Y;    
    data_structure.xdata=XY;
    
    %iterate thru clusters
%     for i= 2168:2169
        if ~isempty(RFprops.error) % ignore clusters with any error message
            return;
        end
        
        %if needed, substract the mean
%         RFi=RF(:,:,:,i);
%         if substract_mean
%             RFi_mean=median(RFi,3);  
%             RFi_nomean = RFi-RFi_mean;
%         else
%             RFi_nomean = RFi;
%         end           
        RFi_nomean = rf;
        
        model=Gaussian_Sum_Rot_center_surround_sep([RFprops.center,RFprops.surround],XY);
        
        %make the templates of the oval corressponding to positive and
        %negative gaussians             
        %center, or the stronger component       
        if RFprops.center_exist==1
            center_mask = make_ellipse_mask(nrow,ncol, RFprops.center(2:end));
            
            A=RFprops.center(1);
            mxc=RFprops.center(2);
            myc=RFprops.center(3);
            sigma1=RFprops.center(4);
            sigma2=RFprops.center(5);
            center_thresh = A*exp(   -((sigma1-mxc).^2/(2*sigma1^2) + (sigma2-myc).^2/(2*sigma2^2) )    );
            
            if A>0
                center_vals=model.*(model>0);
            else
                center_vals=abs(model.*(model<0));
            end
            sumcentr=sum(sum(center_vals));
            xcw=sum(sum(Xt.*center_vals))/sumcentr;
            ycw=sum(sum(Yt.*center_vals))/sumcentr;    
%             [mincenter,mincenter_t]=min(RFprops.center_dynamics*RFprops.center(1));
%             [maxcenter,maxcenter_t]=max(RFprops.center_dynamics*RFprops.center(1));
%             time_max_min=abs(mincenter_t-maxcenter_t);
        end

        %surround, or the weaker component
        if RFprops.surround_exist==1
            surround_mask = make_ellipse_mask(nrow,ncol, RFprops.surround(2:end));
            surround_mask=surround_mask-center_mask;
            surround_mask(surround_mask<0)=0;
            
            A=RFprops.surround(1);
            mxs=RFprops.surround(2);
            mys=RFprops.surround(3);
            sigma1=RFprops.surround(4);
            sigma2=RFprops.surround(5);
            surround_thresh = A*exp(   -((sigma1-mxs).^2/(2*sigma1^2) + (sigma2-mys).^2/(2*sigma2^2) )    );
            
            if A>0
                surr_vals=model.*(model>0);
            else
                surr_vals=abs(model.*(model<0));
            end       
            sumsurr=sum(sum(surr_vals));
            xsw=sum(sum(Xt.*surr_vals))/sumsurr;
            ysw=sum(sum(Yt.*surr_vals))/sumsurr;
            
%             [minsurr,minsurr_t]=min(RFprops.surround_dynamics*RFprops.surround(1));
%             [maxsurr,maxsurr_t]=max(RFprops.surround_dynamics*RFprops.surround(1));  
%             time_max_min=abs(minsurr_t-maxsurr_t);
        end
        
        %if both components exist compute inter-component properties
        if RFprops.center_exist==1 && RFprops.surround_exist==1
            %distance between gaussian centers
            d=sqrt((mxs-mxc)^2+(mys-myc)^2);
            d_weightedcenters=sqrt((xsw-xcw)^2+(ysw-ycw)^2);
            %angle between gaussian centers and horizontal
%             alpha=acos(dot([mxs-mxc, mys-myc],[1,0])/d);
            alpha = atan2(mys-myc, mxs-mxc);
%             alpha_weightedcenters=acos(dot([xsw-xcw, ysw-ycw],[1,0])/d_weightedcenters);
            alpha_weightedcenters = atan2(ysw-ycw, xsw-xcw);
            %time between the time peaks as fitted from gaussians
%             if maxcenter>maxsurr, maxt=maxcenter_t;
%             else, maxt=maxsurr_t;
%             end
%             
%             if mincenter<minsurr, mint=mincenter_t;
%             else, mint=minsurr_t;
%             end
%             dt=abs(maxt-mint);
        else
            d=nan; d_weightedcenters=nan;
            alpha=nan; alpha_weightedcenters=nan;            
        end
       
       data_structure.param = [RFprops.center, RFprops.surround];
    
       if RFprops.center_exist==-1
            data_structure.param(1:6)=0;
       end
       if RFprops.surround_exist==-1
            data_structure.param(7:12)=0;
       end
       

        RFprops.center_mask=center_mask;
        RFprops.surround_mask=surround_mask;                             

        
        RFprops.distance_btw_gaussian_centers=d;
        RFprops.distance_btw_centers_of_mass=d_weightedcenters;
        RFprops.orientation_gaussians=alpha; 
        RFprops.orientation_centers_of_mass=alpha_weightedcenters; 
        

        
        %model without scale params
        RFprops.model=model;          
%     end  
%     propsfile=fullfile(res_folder,'RF_props.mat');
%     save(propsfile,'RFprops');


end

function w2=half_width_from_dynamics(ydyn)
    [ypeak,xpeak]=max(abs(ydyn));
    h_2=ypeak/2;
    h_2_signed=ydyn(xpeak)/2;
    %find the point on the left of peak with y(x)=h_2
    xbefore=find(abs(ydyn(1:xpeak))<h_2,1,'last');
    ybefore=ydyn(xbefore);
    yafter=ydyn(xbefore+1);
    x1=interp1([ybefore, yafter],[xbefore,xbefore+1],h_2_signed,'linear');
    
    %find the point on the right of peak with y(x)=h_2
    xafter=find(abs(ydyn(xpeak:end))<h_2,1,'first');
    xafter=xafter+xpeak-1;
    yafter=ydyn(xafter);
    ybefore=ydyn(xafter-1);
    x2=interp1([ybefore, yafter],[xafter-1,xafter],h_2_signed,'linear');
    w2=x2-x1;
end

function RFprops = set_void_components(RFprops)
    ncl=length(RFprops);
    for i=1:ncl
        %skip components with error message
        if ~isempty(RFprops.error)
            RFprops.center_exist=-1;
            RFprops.surround_exist=-1;
        else
            if ~isfield(RFprops, "center")
                RFprops.center_exist=-1;
            else
                RFprops.center_exist=1;
            end
            if ~isfield(RFprops, "surround")
                RFprops.surround_exist=-1;
            else
                RFprops.surround_exist=1;
            end            
        end
    end
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
    
    mystrel=strel('disk', max(rp1,rp2));
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
        RFprops.id=id_depth(i,1);
        RFprops.depth=id_depth(i,2);
    end    
end




function [xn,yn] = rotate_pts(pt,ptcenter, alpha)
    pttemp=pt-ptcenter;
    xn = pttemp(1)*cos(alpha) - pttemp(2)*sin(alpha);
    yn = pttemp(1)*sin(alpha) + pttemp(2)*cos(alpha);
    xn=xn+ptcenter(1);
    yn=yn+ptcenter(2);
end
