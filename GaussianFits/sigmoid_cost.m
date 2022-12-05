function y=sigmoid_cost(x,lowbd,highbd)
    mincost=0;
    maxcost=1000;
    steepness=10;
    midbd=lowbd+0.5*(highbd-lowbd);    
    y= maxcost./(1+exp(-steepness*(x-midbd)))+mincost;
end