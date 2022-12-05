%% use this script to plot properties of RF
% computed from fitting two gaussians to RF
% (first run RF_gaussian_mixture.m and gaussian_mixture_RF_properties.m)
% Structure RFprops should be loaded 

% O.Symonova
% 11.02.2021


%% folder to store figure
res_folder='C:\DATA\Data_for_Olga\retina\';
substract_mean=1;
%% if depth info is available provide depth file
%     depth_file = 'C:\DATA\EPhys\silicone_probe\291020\cluster_info.tsv';
    depth_file=[];
    if ~isempty(depth_file)
        with_depth=1;
    else
        with_depth=0;
    end
    
    ncl=length(RFprops);
    noerrcl=0;
    [nrow,ncol]=size(RFprops(1).model);
    
%% count how many clusters with at least one component: this is the nuber of rows in the table
    for i=1:ncl
        if RFprops(i).center_exist==1 || RFprops(i).surround_exist==1
            noerrcl=noerrcl+1;
        end
    end
    
    miny=1000; maxy=1; maxw1=0;
    minx=1000; maxx=1; maxw2=0;
    
%% create a data and model image for the three frames around the max variance time

    images_data=zeros(nrow,ncol,ncl);
    images_model=zeros(nrow,ncol,ncl);
    %create the grid to compute model values
    [X, Y] = meshgrid(1:ncol,1:nrow);
    XY=[];
    XY(:,:,1) = X;
    XY(:,:,2) = Y;    
    data_structure.xdata=XY;
    
    for i=1:ncl
        if ~isempty(RFprops(i).error)
            continue;
        end
        %substract meadian to bring baseline to 0
        RFi=RF(:,:,:,i);
        if substract_mean
            RFi_mean=median(RFi,3);  
            RFi_nomean = RFi-RFi_mean;
        else
            RFi_nomean = RFi;
        end     
        %indices around max varaince
        tid=RFprops(i).maxvar_ind;
        tst=max(1,tid-1);
        tend=min(nhist,tid+1);
        %RF around max varians
        imdata = RFi_nomean(:,:,tst:tend);
        imdata =median(imdata,3);
        images_data(:,:,i)=imdata;
    end
    
    %% get min-max values for cropping
    for i=1:ncl
        if ~isempty(RFprops(i).error)
            continue;
        end
        miny=min([miny, RFprops(i).center(3), RFprops(i).surround(3)]);
        maxy=max([maxy, RFprops(i).center(3), RFprops(i).surround(3)]);
        minx=min([minx, RFprops(i).center(2), RFprops(i).surround(2)]);
        maxx=max([maxx, RFprops(i).center(2), RFprops(i).surround(2)]);
        maxw1=max([maxw1, RFprops(i).center(5), RFprops(i).surround(5)]);
        maxw2=max([maxw2, RFprops(i).center(6), RFprops(i).surround(6)]);
    end
    maxw=max(maxw1,maxw2);
    mincropy = max(1,round(miny-maxw));  maxcropy = min(round(maxy+maxw),nrow);
    mincropx = max(1,round(minx-maxw));  maxcropx = min(round(maxx+maxw),ncol);
    
    minval=realmax;
    maxval=realmin;
    
    maxav=realmin;
    minav=realmax;
    maxcomponent=realmin;
    mincomponent=realmax;
    maxsize=0;
    maxdyn=0;
            
    for i=1:ncl
        if ~isempty(RFprops(i).error)
            continue;
        end
        visfield=RFprops(i).model(mincropy:maxcropy,mincropx:maxcropx);
        minval=min(minval,min(visfield(:)));
        maxval=max(maxval,max(visfield(:))); 
        maxdyn=max([maxdyn,max(abs(RFprops(i).center_dynamics*RFprops(i).center(1))), max(abs(RFprops(i).surround_dynamics*RFprops(i).surround(1)))]);
        maxav=max([maxav,max(abs(RFprops(i).center_av_value_data)), max(abs(RFprops(i).center_av_value_fit))]);
        maxav=max([maxav,max(abs(RFprops(i).surround_av_value_data)), max(abs(RFprops(i).surround_av_value_fit))]);
        maxcomponent=max([maxcomponent,max(abs(RFprops(i).center_same_sign_data)), max(abs(RFprops(i).center_same_sign_fit))]);
        maxcomponent=max([maxcomponent,max(abs(RFprops(i).surround_same_sign_data)), max(abs(RFprops(i).surround_same_sign_fit))]);
        maxsize=max([maxsize,max(RFprops(i).center_same_sign_size_data), max(RFprops(i).center_same_sign_size_fit),...
            max(RFprops(i).surround_same_sign_size_data), max(RFprops(i).surround_same_sign_size_fit)]);
    end
    
    %if depth is provided sort clusters by depth
    if with_depth
        ids_depth = get_clusters_depth(depth_file);
        id_depth=[1:ncl; ids_depth(:,2)']';
        id_depth=sortrows(id_depth,2);
    end
    
    bluredmap=bluered();
   
    plot1_hight=10; %height of each subplot in some units - later will be mulitiplied by 10 to convert to pixels
    plot1_width=10;
    marg=plot1_hight/5;
    labelw=plot1_hight/2;
    figurew=labelw+8*plot1_width+8*marg;%labelw; %width of the figure
    figureh=plot1_hight*noerrcl+2*marg+labelw; %height of the figure

   
    f= figure;
    set(f,'PaperPosition',[0,0,figurew, figureh]); %size of the figure to save is larger than than screen
    set(f,'Visible','off');    
    
    %values for margins and spaces btw subplots
    l1=0;
    pw=0.13;
    tw=pw/2;
    mm=pw/5;
    l2=tw;
    l3=l2+pw+mm;
    l4=l3+pw+mm;
    l5=l4+pw+mm;
    l6=l5+pw+mm;
    l7=l6+pw+mm;
    
    heightit=1/noerrcl; 
    mv=heightit/6;
    height_p=1-2*mv;
    heighti=height_p/noerrcl - mv; 
    bottomi=1-mv-heighti;
    firstcluster=1;
    
    if with_depth
        indices=id_depth(:,1);
    else
        indices=1:ncl;
    end
    
    for i_i=1:ncl
        i=indices(i_i);        
        
        if RFprops(i).center_exist==-1 && RFprops(i).surround_exist==-1
            continue;
        end                
              
        %depth label   
        subplot('Position',[l1,bottomi,tw,heighti]); 
        axis off;
        [max_center,peak_time]=max(RFprops(i).center_dynamics);
        peaktime_st=max(1,peak_time-10);
        peaktime_en=min(length(RFprops(i).center_dynamics),peak_time-10);
        max_center=abs(max_center*RFprops(i).center(1));
        max_surround=abs(RFprops(i).surround_dynamics(peaktime_st:peaktime_en)*RFprops(i).surround(1));
        if with_depth            
            %ii = id_depth(i,1);
            depthi=id_depth(i,2);
            text(0.1,0.75,{[num2str(depthi,'%d'),'(',num2str(i,'%d'),')'],...
            ['max/mean var: ',num2str(RFprops(i).max_maean_variance,'%.1f')],...
            ['s-c correlation: ',num2str(RFprops(i).surround_center_corr,'%.2f')],...
            ['s-c strength: ',num2str(max_surround/max_center,'%.4f')]}, 'FontSize',50);
            %text(0.1,0.5,[num2str(depthi,'%d'),'(',num2str(i,'%d'),')'], 'FontSize',50);
        else
             %text(0.1,0.5,['cl ',num2str(i,'%d')], 'FontSize',50);
             text(0.1,0.75,{['cl ',num2str(i,'%d')],...
             ['max/mean var: ',num2str(RFprops(i).max_maean_variance,'%.1f')],...
             ['s-c correlation: ',num2str(RFprops(i).surround_center_corr,'%.2f')],...
             ['s-c strength: ',num2str(max_surround/max_center,'%.4f')]}, 'FontSize',50);
        end
        
        
        %max var data
        ax = subplot('Position',[l2, bottomi, pw, heighti]);
                
        imdata=images_data(:,:,i);
        %set [-1,1]
        maxval=max(max(imdata(:)),abs(min(imdata(:))));
%         imdata=imdata/maxval;
        imagesc(imdata(mincropy:maxcropy, mincropx:maxcropx),[-maxval,maxval]);
        colormap(bluredmap);
        hold on;
        if RFprops(i).center_exist==1
            ell_pos=RFprops(i).center;             
            plot_ellipse(ell_pos(2)-mincropx+1,ell_pos(3)-mincropy+1,ell_pos(4)*2,ell_pos(5)*2,-ell_pos(6),[1,0,0],0.0);
            xlim([mincropx,maxcropx]);
            ylim([mincropy,maxcropy]);
            hold on;
        end
        if RFprops(i).surround_exist==1
            ell_neg=RFprops(i).surround;
            plot_ellipse(ell_neg(2)-mincropx+1,ell_neg(3)-mincropy+1,ell_neg(4)*2,ell_neg(5)*2,-ell_neg(6),[0,0,1],0.0);
            xlim([mincropx,maxcropx]);
            ylim([mincropy,maxcropy]);
            hold off;
        end
%         xlim([mincropx,maxcropx]);
%         ylim([mincropy,maxcropy]);
        axis equal;
        axis tight;
        xticks([]); yticks([]); 
        
%         colorbar;
        if firstcluster
            ax2=ax;
        end
                
%         %max model fit
%         ax = subplot('Position',[l2, bottomi, pw, heighti]);
% %         imagesc(RFprops(i).model_max(mincropy:maxcropy, mincropx:maxcropx),[minval,maxval]);
%         imagesc(RFprops(i).model_max(mincropy:maxcropy, mincropx:maxcropx));
%         axis equal;
%         axis tight;
%         xticks([]); yticks([]);        
%         colormap(bluredmap);
% %         colorbar;
%         if firstcluster
%             ax2=ax;
%         end
        

        %max var model
        ax = subplot('Position',[l3, bottomi, pw, heighti]);
                
        imdata=RFprops(i).model;
        %set [-1,1]
        maxval=max(abs(imdata(:)));
%         imdata=imdata/maxval;
        imagesc(imdata(mincropy:maxcropy, mincropx:maxcropx),[-maxval,maxval]);
        colormap(bluredmap);
        axis equal;
        axis tight;
        xticks([]); yticks([]);        
                
        
        %component dynamics
        ax = subplot('Position',[l4, bottomi, pw, heighti]);
        plot(RFprops(i).center_dynamics*RFprops(i).center(1),'r'), hold on, 
        plot(RFprops(i).surround_dynamics*RFprops(i).surround(1),'b');
%         ylim([-maxdyn,maxdyn]); %normalize the plot to global max and min
        xticks([]); 
        maxval=max(abs([RFprops(i).center_dynamics'*RFprops(i).center(1),RFprops(i).surround_dynamics'*RFprops(i).surround(1)]));
        ylim([-maxval-maxval*0.1,maxval+maxval*0.1]);
        yticks([-maxval,0,maxval]); 
        yticklabels({});
        text(ax, [-0.5,-0.5,-0.5], [-maxval,0,maxval],{num2str(-maxval,'%.1f'),'0',num2str(maxval,'%.1f')},'FontSize',32,'HorizontalAlignment','right');      
        if firstcluster
            ax3=ax;
        end
        
        %average component values
        ax = subplot('Position',[l5, bottomi, pw, heighti]);
        if RFprops(i).center_exist==1, plot(RFprops(i).center_av_value_data,'r'), hold on, plot(RFprops(i).center_av_value_fit,'--r'), hold on; end
        if RFprops(i).surround_exist==1, plot(RFprops(i).surround_av_value_data,'b'), hold on, plot(RFprops(i).surround_av_value_fit,'--b'); end
%         ylim([-maxav, maxav]); 
        xticks([]);         
        maxval=max([max(abs(RFprops(i).center_av_value_data)), max(abs(RFprops(i).center_av_value_fit)),...
            max(abs(RFprops(i).surround_av_value_data)), max(abs(RFprops(i).surround_av_value_fit))]);
        ylim([-maxval-maxval*0.1,maxval+maxval*0.1]);
        yticks([-maxval,0,maxval]); 
        yticklabels({});
        text(ax, [-0.5,-0.5,-0.5], [-maxval,0,maxval],{num2str(-maxval,'%.1f'),'0',num2str(maxval,'%.1f')},'FontSize',32,'HorizontalAlignment','right');      
        if firstcluster
            ax4=ax;
        end
        
        ax = subplot('Position',[l6, bottomi, pw, heighti]);
        if RFprops(i).center_exist==1, plot(RFprops(i).center_same_sign_data,'r'), hold on, plot(RFprops(i).center_same_sign_fit,'--r'), hold on; end
        if RFprops(i).surround_exist==1, plot(RFprops(i).surround_same_sign_data,'b'), hold on, plot(RFprops(i).surround_same_sign_fit,'--b'), hold on; end
%         ylim([-maxcomponent,maxcomponent]); 
        xticks([]); %yticks([]); 
        maxval=max([max(abs(RFprops(i).center_same_sign_data)), max(abs(RFprops(i).center_same_sign_fit)),...
            max(abs(RFprops(i).surround_same_sign_data)), max(abs(RFprops(i).surround_same_sign_fit))]);
        ylim([-maxval-maxval*0.1,maxval+maxval*0.1]);
        yticks([-maxval,0,maxval]); 
        yticklabels({});
        text(ax, [-0.5,-0.5,-0.5], [-maxval,0,maxval],{num2str(-maxval,'%.1f'),'0',num2str(maxval,'%.1f')},'FontSize',32,'HorizontalAlignment','right');      
        if firstcluster
            ax5=ax;
        end
        
        
        ax = subplot('Position',[l7, bottomi, pw, heighti]);
        if RFprops(i).center_exist==1, plot(RFprops(i).center_same_sign_size_data,'r'), hold on, plot(RFprops(i).center_same_sign_size_fit,'--r'), hold on; end
        if RFprops(i).surround_exist==1, plot(RFprops(i).surround_same_sign_size_data,'b'), hold on, plot(RFprops(i).surround_same_sign_size_fit,'--b'); end
        %ignore the limit for now, noisy clusters dominate the size
%         ylim([0,maxsize]); 
        maxi=max([RFprops(i).center_same_sign_size_data, RFprops(i).center_same_sign_size_fit, RFprops(i).surround_same_sign_size_data, RFprops(i).surround_same_sign_size_fit]);
        xticks([]);
        ylim([-1,maxi]);
        yticks([0,maxi]); 
        yticklabels({});        
        text(ax,[-0.5,-0.5],[0,maxi],{'0',num2str(maxi)},'FontSize',32,'HorizontalAlignment','right');         
%         ytickformat('%d');
%         ax.YAxisLocation = 'right';
        if firstcluster
            ax6=ax;
        end
        
        firstcluster=0;
        bottomi=bottomi-mv-heighti;

    end

title(ax2, 'amplitutde','FontSize',50);
title(ax3, 'model','FontSize',50);
title(ax4, 'average total','FontSize',50);
title(ax5, 'average same sign','FontSize',50);
title(ax6, 'size','FontSize',50);

pngfile=fullfile(res_folder,'gaussian_properties_new.png');
print(gcf,'-r20','-dpng',pngfile);

        
            
          
        
%     end
    
    
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
    
