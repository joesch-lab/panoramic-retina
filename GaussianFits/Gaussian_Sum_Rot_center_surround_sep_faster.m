function F = Gaussian_Sum_Rot_center_surround_sep_faster(params,datamix)
%% data_structure.param = [Amp1, x1, y1, wx1, wy1, fi1, Amp2, wx2, wy2, fi2]

    g1=[params(1),datamix.g1];
    g2=params(2:end);
    xdata=datamix.data(:,2:end);
    g1space = datamix.data(:,1);
    X2=xdata(:,1)-g2(2);
    Y2=xdata(:,2)-g2(3);    
    XY2=[X2,Y2];
    
    F1 = g1(1)*g1space;    
    
    %2nd gaussian
    rotM2=rotation_matrix(g2(6));   
    xdatarot2=XY2*rotM2;    
    xdatarot2(:,1)=xdatarot2(:,1)+g2(2);
    xdatarot2(:,2)=xdatarot2(:,2)+g2(3);
    F2 = g2(1)*exp(   -((xdatarot2(:,1)-g2(2)).^2/(2*g2(4)^2) + (xdatarot2(:,2)-g2(3)).^2/(2*g2(5)^2) )    );
        
    F=F1+F2;
end