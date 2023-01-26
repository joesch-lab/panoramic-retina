
function cents = polar_to_visual(polar, alpha0, theta0)
lambda_all = polar(:,1);
phi_all = polar(:,2);

lambda = wrapTo2Pi(lambda_all+pi/2);
[phi,lambda] = invert_sphere(phi_all,lambda);
[x,y,z] = sphere2cart(phi,lambda);
[x,y,z] = rotate_axis(x,y,z, theta0, alpha0);
[phi,lambda] = cart2sphere(x,y,z);
[x,y] = sinusoidal_transform(phi,lambda);
cents = [x y];

    function [x,y] = sinusoidal_transform(phi, lambda)
        lambda0 = 0; %mid longitude of plot (GMT or IDL is in the center?)
        
        x = (lambda-lambda0).*cos(phi);
        y = phi;
    end

    function [phi, lambda] = invert_sphere(phi,lambda)
        phi = -phi;
        lambda = wrapToPi(lambda + pi);
    end

    function [x,y,z] = sphere2cart(phi,lambda)
        x = -cos(phi).*sin(lambda);
        z = sin(phi);
        y = cos(phi).*cos(lambda);
    end

    function [phi,lambda] = cart2sphere(x,y,z)
        phi = atan(z./sqrt(x.^2+y.^2));
        lambda = atan2(-x,y);
    end

    %theta0 - final azimuth of optical axis (+- 60)
    %alpha0 - final elevation of optical axis ( +22)
    function [x2,y2,z2] = rotate_axis(x0,y0,z0,theta0, alpha0)
        %initially optic axis is vertical (+z)
        %first rotate -90+alpha0 about x axis
        % so that axis is at correct elevation
        rot = -pi/2+alpha0;
        x1 = x0;
        y1 = y0*cos(rot) - z0*sin(rot);
        z1 = z0*cos(rot) + y0*sin(rot);
        
        %then rotate about z axis such that azimuth ends up at theta0
        rot = -theta0;
        z2 = z1;
        x2 = x1*cos(rot) - y1*sin(rot);
        y2 = y1*cos(rot) + x1*sin(rot);
    end
end

