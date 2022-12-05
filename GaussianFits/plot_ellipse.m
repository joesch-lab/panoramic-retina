function plot_ellipse(xc, yc, wx, wy, angle, color, transparency)
    %a: width in pixels
    %b: height in pixels
    %xc: horizontal center
    %yc: vertical center
    %angle: orientation ellipse in radians
    %color: color code (e.g., 'r' or [0.4 0.5 0.1])
%     a=std_laten(1)*pca_axis(1,:);
%     b=std_laten(2)*pca_axis(2,:);    
    
    
%     b = pca_ellipse_axis(2,:);
%     a=pca_ellipse_axis(1,:);
%     alen=sqrt(a(1)^2+a(2)^2);
%     blen=sqrt(b(1)^2+b(2)^2);
%     angle = -acos(dot(a,[1,0])/(alen));    
    
    a=wx*[1,0];
    b=wy*[0,1];
    alen=wx;
    blen=wy;

    r=0:0.1:2*pi+0.1;
    p=[(alen*cos(r))' (blen*sin(r))'];
    alpha=[cos(angle) -sin(angle)
           sin(angle) cos(angle)];

     p1=p*alpha;
    patch(xc+p1(:,1),yc+p1(:,2),color,'FaceAlpha',transparency,'EdgeColor',color, 'LineWidth', 2);
end
