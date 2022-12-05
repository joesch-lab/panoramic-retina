function RFprops = gaussian_mixture_RF_properties_circular(RF,RFprops)
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
    [nrow,ncol, nhist, ncl]=size(RF);
    
    %create the grid to compute model values
    [X, Y] = meshgrid(1:ncol,1:nrow);
    Xt = X;
    Yt = Y;
    XY=[];
    XY(:,:,1) = X;
    XY(:,:,2) = Y;    
    data_structure.xdata=XY;

    %make a circular mask    
    rf_center = floor([ncol,nrow]/2)+1;
    maxR=floor(min([ncol,nrow])/2);
    circ_mask = ((Xt-rf_center(1)).^2+(Yt-rf_center(2)).^2).^0.5<maxR;
    
    %iterate thru clusters
    for i=1:ncl
        if mod(i,500)==0
            disp(['Processing cluser ', num2str(i), ' out of ', num2str(ncl),'.'])
        end
        if ~isempty(RFprops(i).error) % ignore clusters with any error message
            continue;
        end
        
        %if needed, substract the mean
        RFi=RF(:,:,:,i);
        if substract_mean
            RFi_mean=median(RFi,3);  
            RFi_nomean = RFi-RFi_mean;
        else
            RFi_nomean = RFi;
        end                
        
        if RFprops(i).center_exist==1 && RFprops(i).surround_exist==1
            model=Gaussian_Sum_Rot_center_surround_sep([RFprops(i).center,RFprops(i).surround],XY);
        elseif RFprops(i).center_exist==1 
             model=Gaussian_Sum_Rot_center(RFprops(i).center,XY);
        else
            continue;
        end
        
        %make the templates of the oval corressponding to positive and
        %negative gaussians             
        %center, or the stronger component       
        if RFprops(i).center_exist==1
            center_mask = make_ellipse_mask(nrow,ncol, RFprops(i).center(2:end));
            
            A=RFprops(i).center(1);
            mxc=RFprops(i).center(2);
            myc=RFprops(i).center(3);
            sigma1=RFprops(i).center(4);
            sigma2=RFprops(i).center(5);
            center_thresh = A*exp(   -((sigma1-mxc).^2/(2*sigma1^2) + (sigma2-myc).^2/(2*sigma2^2) )    );
            
            if A>0
                center_vals=model.*(model>0);
            else
                center_vals=abs(model.*(model<0));
            end
            sumcentr=sum(sum(center_vals));
            xcw=sum(sum(Xt.*center_vals))/sumcentr;
            ycw=sum(sum(Yt.*center_vals))/sumcentr;    
            [mincenter,mincenter_t]=min(RFprops(i).center_dynamics*RFprops(i).center(1));
            [maxcenter,maxcenter_t]=max(RFprops(i).center_dynamics*RFprops(i).center(1));
            time_max_min=abs(mincenter_t-maxcenter_t);
        end

        %surround, or the weaker component
        if RFprops(i).surround_exist==1
            surround_mask = make_ellipse_mask(nrow,ncol, RFprops(i).surround(2:end));
            surround_mask=surround_mask-center_mask;
            surround_mask(surround_mask<0)=0;
            
            A=RFprops(i).surround(1);
            mxs=RFprops(i).surround(2);
            mys=RFprops(i).surround(3);
            sigma1=RFprops(i).surround(4);
            sigma2=RFprops(i).surround(5);
            surround_thresh = A*exp(   -((sigma1-mxs).^2/(2*sigma1^2) + (sigma2-mys).^2/(2*sigma2^2) )    );
            
            if A>0
                surr_vals=model.*(model>0);
            else
                surr_vals=abs(model.*(model<0));
            end       
            sumsurr=sum(sum(surr_vals));
            xsw=sum(sum(Xt.*surr_vals))/sumsurr;
            ysw=sum(sum(Yt.*surr_vals))/sumsurr;
            
            [minsurr,minsurr_t]=min(RFprops(i).surround_dynamics*RFprops(i).surround(1));
            [maxsurr,maxsurr_t]=max(RFprops(i).surround_dynamics*RFprops(i).surround(1));  
            time_max_min=abs(minsurr_t-maxsurr_t);
        end
        
        %if both components exist compute inter-component properties
        if RFprops(i).center_exist==1 && RFprops(i).surround_exist==1
            %distance between gaussian centers
            d=sqrt((mxs-mxc)^2+(mys-myc)^2);
            d_weightedcenters=sqrt((xsw-xcw)^2+(ysw-ycw)^2);
            %angle between gaussian centers and horizontal
%             alpha=acos(dot([mxs-mxc, mys-myc],[1,0])/d);
            alpha = atan2(mys-myc, mxs-mxc);
%             alpha_weightedcenters=acos(dot([xsw-xcw, ysw-ycw],[1,0])/d_weightedcenters);
            alpha_weightedcenters = atan2(ysw-ycw, xsw-xcw);
            %time between the time peaks as fitted from gaussians
            if maxcenter>maxsurr, maxt=maxcenter_t;
            else, maxt=maxsurr_t;
            end
            
            if mincenter<minsurr, mint=mincenter_t;
            else, mint=minsurr_t;
            end
            dt=abs(maxt-mint);

            %% new center-of-mass location and orientation on circle-cropped RF
            % it is a weighted sum of pixels of the model
            sumcentr=sum(center_vals(circ_mask));
            Xw=Xt.*center_vals;
            xcw_circ=sum(Xw(circ_mask))/sumcentr;
            Yw=Yt.*center_vals;
            ycw_circ=sum(Yw(circ_mask))/sumcentr;
        
          
            sumsurr=sum(surr_vals(circ_mask));
            Xw=Xt.*surr_vals;
            xsw_circ=sum(Xw(circ_mask))/sumsurr;
            Yw=Yt.*surr_vals;
            ysw_circ=sum(Yw(circ_mask))/sumsurr;    
            
            alpha_weightedcenters_circ = atan2(ysw_circ-ycw_circ, xsw_circ-xcw_circ);
    
            %% new center-of-mass location and orientation on mask-cropped RF
            % it is a mean of data RF within center and surround masks
            center_data_weighted = RFprops(i).RFaroundpeak.*center_mask;
            center_data_weighted(sign(center_data_weighted)~=sign(RFprops(i).center(1))) =0;
            center_data_weighted=abs(center_data_weighted);
            sumcentr=sum(center_data_weighted(:));
            Xw=Xt.*center_data_weighted;
            xcw_data_mask=sum(Xw(:))/sumcentr;
            Yw=Yt.*center_data_weighted;
            ycw_data_mask=sum(Yw(:))/sumcentr;
    
            surround_data_weighted = RFprops(i).RFaroundpeak.*surround_mask;
            surround_data_weighted(sign(surround_data_weighted)~=sign(RFprops(i).surround(1))) =0;
            surround_data_weighted=abs(surround_data_weighted);
            sumsurr=sum(surround_data_weighted(:));
            Xw=Xt.*surround_data_weighted;
            xsw_data_mask=sum(Xw(:))/sumsurr;
            Yw=Yt.*surround_data_weighted;
            ysw_data_mask=sum(Yw(:))/sumsurr;
            
            alpha_data_mask = atan2(ysw_data_mask-ycw_data_mask, xsw_data_mask-xcw_data_mask);
        else
            d=nan; d_weightedcenters=nan;
            alpha=nan; alpha_weightedcenters=nan;  
            alpha_weightedcenters_circ=nan; alpha_data_mask=nan;            
        end

        %% save orientation of the min-max line computer on the RF
        [RFijmax,tmaxind] = max(RFi_nomean(:));
        [RFijmin,tminind] = min(RFi_nomean(:));
        [ypeakp, xpeakp, tpeakp] = ind2sub(size(RFi_nomean),tmaxind);            
        [ypeakn, xpeakn, tpeakn] = ind2sub(size(RFi_nomean),tminind); 
        if RFprops(i).center(1)>0 %on-cell
            rot_c_s_line = atan2(ypeakn-ypeakp, xpeakn-xpeakp);
        else
            rot_c_s_line = atan2(ypeakp-ypeakn, xpeakp-xpeakn);
        end
           
       data_structure.param = [RFprops(i).center, RFprops(i).surround];
    
       if RFprops(i).center_exist==-1
            data_structure.param(1:6)=0;
       end
       if RFprops(i).surround_exist==-1
            data_structure.param(7:12)=0;
       end

       %% compute the R2 values        
       R2_c=nan;
       R2_cs=nan;
       dR2=nan;

       if RFprops(i).center_exist==1
            %compute the fit of the model if only the center was fit
            model_center_only = Gaussian_Rot_center(RFprops(i).center,XY); 
            %values of data not explained by the model
            imdiff_c=model_center_only - RFprops(i).RFaroundpeak;
            %bounding box of the center mask to compute R2 value in the box
            
            [c_row,c_col]=find(center_mask>0);
            mincol=min(c_col);
            maxcol=max(c_col);
            minrow=min(c_row);
            maxrow=max(c_row);
    
            if RFprops(i).surround_exist               
                %get the image of the model center-surround
                model_center_surround = Gaussian_Sum_Rot_center_surround_sep([RFprops(i).center, RFprops(i).surround],XY);   
                %data values not explained by the model
                imdiff_cs=model_center_surround-RFprops(i).RFaroundpeak;
                
                %bounding box of the center-surround masks to compute R2 
                %compute measure of fit on the area cropped to the RF        
                [s_row,s_col]=find(surround_mask>0);
                mincol=min(mincol,min(s_col));
                maxcol=max(maxcol,max(s_col));
                minrow=min(minrow,min(s_row));
                maxrow=max(maxrow,max(s_row));
                
                %crop model and data to the RF
                imdiff_cs=imdiff_cs(minrow:maxrow,mincol:maxcol);
                imdiff_c=imdiff_c(minrow:maxrow,mincol:maxcol);
                data_cropped = RFprops(i).RFaroundpeak(minrow:maxrow,mincol:maxcol); 
                
                ymean=mean(data_cropped(:));
                SStot=sum((data_cropped(:)-ymean).^2);
                diffcol=imdiff_cs(:);        
                SSres_cs=(sum(diffcol.^2));
                R2_cs=1-SSres_cs/SStot;
                
                diffcol=imdiff_c(:);
                SSres_c=(sum(diffcol.^2));
                R2_c=1-SSres_c/SStot;
                
                dR2=R2_cs-R2_c;    
            else
                %crop model and data to the RF        
                imdiff_c=imdiff_c(minrow:maxrow,mincol:maxcol);
                data_cropped = RFprops(i).RFaroundpeak(minrow:maxrow,mincol:maxcol); 
        
                ymean=mean(data_cropped(:));
                SStot=sum((data_cropped(:)-ymean).^2);        
                
        
                diffcol=imdiff_c(:);
                SSres_c=(sum(diffcol.^2));
                R2_c=1-SSres_c/SStot;
                
                R2_cs=0;
                dR2=R2_cs-R2_c;
            end
        end
        RFprops(i).R2_c=R2_c;
        RFprops(i).R2_cs=R2_cs;
        RFprops(i).dR2=dR2;

           
        %% the values below are computed from the masks of gaussians,
               
        %max pos and negative amplitudes in data and model
        c_max_value_fit=nan(1,nhist);
        c_max_value_data=nan(1,nhist);
        s_max_value_fit=nan(1,nhist);
        s_max_value_data=nan(1,nhist);
       
        %average value of all pixels in pos component in data and model
        c_av_value_fit=nan(1,nhist);
        c_av_value_data=nan(1,nhist);
        
        %average value of all pixels in neg component in data and model
        s_av_value_fit=nan(1,nhist);
        s_av_value_data=nan(1,nhist);
        
        %average value of only positive pixels in pos component in data and model
        c_same_sign_fit=nan(1,nhist);
        c_same_sign_data=nan(1,nhist);
        
        %average value of only neg pixels in neg component in data and model
        s_same_sign_fit=nan(1,nhist);
        s_same_sign_data=nan(1,nhist);        
        
        %size (number of pos pixels) in pos component in data and model
        c_same_sign_size_fit=nan(1,nhist);
        c_same_sign_size_data=nan(1,nhist);
        
        %size (number of neg pixels) in ned component in data and model                
        s_same_sign_size_fit=nan(1,nhist);
        s_same_sign_size_data=nan(1,nhist);
        maxvari=0;
        maxvari_ind=0;
        
        for hi=1:nhist
            %data
            datai=RFi_nomean(:,:,hi);                     
            
            %find in which time sample there is the highest variance in RF
            vari=var(datai(:));
            if vari>maxvari
                maxvari=vari;
                maxvari_ind=hi;
            end
            
            %compute the model
            weightx=[RFprops(i).center_dynamics(hi), RFprops(i).surround_dynamics(hi)];
            F = Gaussian_Mixture_Rot_with_weights(weightx,data_structure);
                        
            
            %1st component
            if RFprops(i).center_exist==1
                 %function value at 1std from peak
                 scale_c=RFprops(i).center_dynamics(hi);
                 %fpos_thresh=scale_pos*pos_thresh;
                 fc_thresh=0;
                center_sign = sign(RFprops(i).center(1)*scale_c);
                %get all the values from data
                %data - 1st component
                datadata=datai.*center_mask;
                nnzvals=nonzeros(datadata);
                if ~isempty(nnzvals)
                    same_sign_vals=nnzvals(center_sign*nnzvals>fc_thresh);
    
                    c_max_value_data(hi)=center_sign*max(center_sign*nnzvals);            
                    c_av_value_data(hi)=mean(nnzvals);
                    if ~isempty(same_sign_vals)
                        c_same_sign_data(hi)=mean(same_sign_vals); 
                        c_same_sign_size_data(hi)=length(same_sign_vals);
                    end
                end
    
                %get all the values from the gaussians
                %model - 1st component
                fitvals=F.*center_mask;
                nnzvals=nonzeros(fitvals);
                if ~isempty(nnzvals)
                    same_sign_vals=nnzvals(center_sign*nnzvals>fc_thresh);
    
                    c_max_value_fit(hi)=center_sign*max(center_sign*nnzvals);           
                    c_av_value_fit(hi)=mean(nnzvals);
                    if ~isempty(same_sign_vals)
                        c_same_sign_fit(hi)=mean(same_sign_vals); 
                        c_same_sign_size_fit(hi)=length(same_sign_vals);                                        
                    end
                end
            end
            
            %2nd component
            if RFprops(i).surround_exist==1  
                %function value at 1std from peak
                scale_s=RFprops(i).surround_dynamics(hi);
                %fneg_thresh=scale_neg*neg_thresh;
                fs_thresh=0;
                %sign of the surround
                surround_sign = sign(RFprops(i).surround(1)*scale_s);
                 
                %data - 2nd component
                datadata=datai.*surround_mask;
                nnzvals=nonzeros(datadata);
                if ~isempty(nnzvals)                    
                    same_sign_vals=nnzvals(surround_sign*nnzvals>fs_thresh);
    
                    s_max_value_data(hi)=surround_sign*max(surround_sign*nnzvals);        
                    s_av_value_data(hi)=mean(nnzvals);
                    if ~isempty(same_sign_vals)
                        s_same_sign_data(hi)=mean(same_sign_vals); 
                        s_same_sign_size_data(hi)=length(same_sign_vals);
                    end
                end
    
                %model - negative component
                fitvals=F.*surround_mask;                
                nnzvals=nonzeros(fitvals);
                if ~isempty(nnzvals)  
                    same_sign_vals=nnzvals(surround_sign*nnzvals>fs_thresh);
    
                    s_max_value_fit(hi)=surround_sign*max(surround_sign*nnzvals);          
                    s_av_value_fit(hi)=mean(nnzvals);
                    if ~isempty(same_sign_vals)
                        s_same_sign_fit(hi)=mean(same_sign_vals); 
                        s_same_sign_size_fit(hi)=length(same_sign_vals);
                    end
                end
            end
        end
            
        RFprops(i).maxvar_ind=maxvari_ind;
        RFprops(i).center_mask=center_mask;
        RFprops(i).surround_mask=surround_mask;                             
        
        RFprops(i).center_max_value_data=c_max_value_data;
        RFprops(i).center_max_value_fit=c_max_value_fit;
        RFprops(i).center_av_value_fit=c_av_value_fit;
        RFprops(i).center_av_value_data=c_av_value_data;
        RFprops(i).center_same_sign_fit=c_same_sign_fit;
        RFprops(i).center_same_sign_data=c_same_sign_data;
        RFprops(i).center_same_sign_size_fit=c_same_sign_size_fit;
        RFprops(i).center_same_sign_size_data=c_same_sign_size_data;
        
        RFprops(i).surround_max_value_fit=s_max_value_fit;
        RFprops(i).surround_max_value_data=s_max_value_data;
        RFprops(i).surround_av_value_fit=s_av_value_fit;
        RFprops(i).surround_av_value_data=s_av_value_data;        
        RFprops(i).surround_same_sign_fit=s_same_sign_fit;
        RFprops(i).surround_same_sign_data=s_same_sign_data;        
        RFprops(i).surround_same_sign_size_fit=s_same_sign_size_fit;
        RFprops(i).surround_same_sign_size_data=s_same_sign_size_data;
        
        RFprops(i).distance_btw_gaussian_centers=d;
        RFprops(i).distance_btw_centers_of_mass=d_weightedcenters;
        RFprops(i).orientation_gaussians=alpha; 
        RFprops(i).orientation_centers_of_mass_old=alpha_weightedcenters; 
        RFprops(i).orientation_centers_of_mass=alpha_weightedcenters_circ; 
        RFprops(i).orientation_centers_of_mass_data_masks=alpha_data_mask;         
        RFprops(i).orientation_min_max_raw=rot_c_s_line;
    
    
        RFprops(i).dt_temporal_peaks=dt;  
    
        
        try 
            RFprops(i).half_width_center = half_width_from_dynamics(RFprops(i).center_dynamics*RFprops(i).center(1));
        catch
        end
        
        %model without scale params
        RFprops(i).model=model;          
    end  
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
        if ~isempty(RFprops(i).error)
            RFprops(i).center_exist=-1;
            RFprops(i).surround_exist=-1;
        else
            if ~any(RFprops(i).center_dynamics)
                RFprops(i).center_exist=-1;
            else
                RFprops(i).center_exist=1;
            end
            if ~any(RFprops(i).surround_dynamics)
                RFprops(i).surround_exist=-1;
            else
                RFprops(i).surround_exist=1;
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
    
    mystrel=strel('disk', max(rp1,rp2),8);
    immystrel=mystrel.Neighborhood;
    immystrel=imresize(immystrel,[2*rp2, 2*rp1],'nearest');
    immystrel=imrotate(immystrel,-alpha_p,'nearest');
    elmask=zeros(nrow,ncol);
    elmask(row_p,col_p)=1;
    elmask=imdilate(elmask,immystrel);
end


