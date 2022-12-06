function F = Gaussian_on_pts(params,xdata)
%% data_structure.param = [Amp1, x1, y1, wx1, wy1, fi1, Amp2, wx2, wy2, fi2]

    g1=params(1:6);    
        
    X=xdata(:,1); X1=X-g1(2); 
    Y=xdata(:,2); Y1=Y-g1(3); 
    XY1=[X1,Y1];    
    
    rotM1=rotation_matrix(g1(6));    
    
    %center
    xdatarot1=XY1*rotM1; 
    xdatarot1(:,1)=xdatarot1(:,1)+g1(2);
    xdatarot1(:,2)=xdatarot1(:,2)+g1(3);
    F = g1(1)*exp(   -((xdatarot1(:,1)-g1(2)).^2/(2*g1(4)^2) + (xdatarot1(:,2)-g1(3)).^2/(2*g1(5)^2) )    );        
end
