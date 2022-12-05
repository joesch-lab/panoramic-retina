function F = Gaussian_Sum_Rot_center_surround_sep(params,xdata)
%% data_structure.param = [Amp1, x1, y1, wx1, wy1, fi1, Amp2, wx2, wy2, fi2]

    g1=params(1:6);
    g2=params(7:12);
    
    [nrow,ncol,nch]=size(xdata);
    X=xdata(:,:,1); X=X(:); X1=X-g1(2); X2=X-g2(2);
    Y=xdata(:,:,2); Y=Y(:); Y1=Y-g1(3); Y2=Y-g2(3);
    XY1=[X1,Y1];
    XY2=[X2,Y2];
    
    rotM1=rotation_matrix(g1(6));    
    
    %1st gaussian        
    xdatarot1=XY1*rotM1; 
    xdatarot1(:,1)=xdatarot1(:,1)+g1(2);
    xdatarot1(:,2)=xdatarot1(:,2)+g1(3);
    F1 = g1(1)*exp(   -((xdatarot1(:,1)-g1(2)).^2/(2*g1(4)^2) + (xdatarot1(:,2)-g1(3)).^2/(2*g1(5)^2) )    );
    
    
    %2nd gaussian
    rotM2=rotation_matrix(g2(6));   
    xdatarot2=XY2*rotM2;    
    xdatarot2(:,1)=xdatarot2(:,1)+g2(2);
    xdatarot2(:,2)=xdatarot2(:,2)+g2(3);
    F2 = g2(1)*exp(   -((xdatarot2(:,1)-g2(2)).^2/(2*g2(4)^2) + (xdatarot2(:,2)-g2(3)).^2/(2*g2(5)^2) )    );
    
    F=F1+F2;
    F=reshape(F,nrow,ncol);
end

function rotM = rotation_matrix(alpha)
    rotM =[cos(alpha), - sin(alpha); sin(alpha),cos(alpha)];
end