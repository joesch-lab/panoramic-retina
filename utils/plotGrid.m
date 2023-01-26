function plotGrid()
lambdas = linspace(-pi,pi,9);
phis = linspace(-pi/2,pi/2,5);
for i = 1:numel(lambdas) %lines of longitude
    lambda = ones(1,100)*lambdas(i);
    phi = linspace(-pi/2,pi/2,100);
    
    [x,y] = sinusoidal_transform(phi,lambda);
    hold on;
    plot(x,y, 'LineWidth', 0.5, 'Color' ,[0.5 0.5 0.5 .5]);
end
for i = 1:numel(phis) %lines of latitude
    phi = ones(1,100)*phis(i);
    lambda = linspace(-pi,pi,100);
    
    [x,y] = sinusoidal_transform(phi,lambda);
    plot(x,y, 'LineWidth', 0.5, 'Color' ,[0.5 0.5 0.5 .5]);
end

    function [x,y] = sinusoidal_transform(phi, lambda)
        lambda0 = 0; %mid longitude of plot (GMT or IDL is in the center?)
        
        x = (lambda-lambda0).*cos(phi);
        y = phi;
    end
end