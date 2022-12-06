function F = GM_oppositesign_scale_params(weightx,data_structure)
%% data_structure.param = [Amp1, x1, y1, wx1, wy1, fi1, Amp2, x2, y2, wx2, wy2, fi2]

    xdataor=data_structure.xdata;
    fixparam=data_structure.param;
    %rotation of the data wrt to the 1st gaussian
    if weightx(1)>0
        fixparam1=fixparam(1:6);
        xdata(:,1)=xdataor(:,1)-fixparam1(2);
        xdata(:,2)=xdataor(:,2)-fixparam1(3);
        xdatarot1(:,1)= xdata(:,1)*cos(fixparam1(6)) - xdata(:,2)*sin(fixparam1(6));
        xdatarot1(:,2)= xdata(:,1)*sin(fixparam1(6)) + xdata(:,2)*cos(fixparam1(6));
        F1 = fixparam1(1)*exp(   -((xdatarot1(:,1)).^2/(2*fixparam1(4)^2) + (xdatarot1(:,2)).^2/(2*fixparam1(5)^2) )    );
    else
        F1=zeros(size(xdataor,1),1);
    end
    %rotation of the data wrt to the 2nd gaussian
    if weightx(2)>0
        fixparam2=fixparam(7:end);
        xdata(:,1)=xdataor(:,1)-fixparam2(2);
        xdata(:,2)=xdataor(:,2)-fixparam2(3);
        xdatarot2(:,1)= xdata(:,1)*cos(fixparam2(6)) - xdata(:,2)*sin(fixparam2(6));
        xdatarot2(:,2)= xdata(:,1)*sin(fixparam2(6)) + xdata(:,2)*cos(fixparam2(6));
        F2 = fixparam2(1)*exp(   -((xdatarot2(:,1)).^2/(2*fixparam2(4)^2) + (xdatarot2(:,2)).^2/(2*fixparam2(5)^2) )    );
    else
        F2=zeros(size(xdataor,1),1);
    end

    F=weightx(1)*F1+(-1)*sign(weightx(1))*weightx(2)*F2;
end

