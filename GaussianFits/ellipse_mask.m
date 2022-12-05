function mask = ellipse_mask(xc, yc, wx, wy, angle, xs, ys)
    %a: width in pixels
    %b: height in pixels
    %xc: horizontal center
    %yc: vertical center
    %angle: orientation ellipse in radians
    % xs, ys from meshgrid

    %remember to use axis xy for plotting

    
mask = zeros(size(xs));

    alen=wx;
    blen=wy;

 
    alpha=[cos(-angle) -sin(-angle)
           sin(-angle) cos(-angle)]; 
       %negative rotation of all the points because we then compare them to the unrotated ellipse
       
       p = [xs(:)-xc ys(:)-yc];
       p1 = p*alpha;
       xs = reshape(p1(:,1), size(mask)) + xc;
       ys = reshape(p1(:,2), size(mask)) + yc;
       
       mask = ((xs-xc).^2/alen^2 + (ys-yc).^2/blen^2) < 1;
end