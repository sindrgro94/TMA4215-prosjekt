% Task 10 - work-precision Diagram: (for Van der Pol equation)
close all
clear all
clc
mu = 10;
f = @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
jac = jacobi('Van der Pol equation',mu);
tInt = [0 mu];
y0 = [2;0];
m = 2; %from y
h0 = 0.1; %can be chosen
options = odeset('RelTol',2.3e-14,'AbsTol',1e-16);
errorEstimator = ode15s(f,tInt,y0,options);
errorEstimator = errorEstimator.y(:,end);%"correct" y at t-end
ourError = zeros(1,10);
ourWork = zeros(1,10);
matlabError15s = zeros(1,10);
matlabWork15s= zeros(1,10);
matlabError23s = zeros(1,10);
matlabWork23s= zeros(1,10);
for tol = 1:10
    [t, y, iflag, nfun, njac] = RKs(f, jac, tInt(1), tInt(2), y0, 10^(-tol), h0,{false});
    stop = find(t==tInt(2));
    ourWork(tol) = nfun+m*njac;
    ourError(tol) = norm(errorEstimator-y(:,stop));
    options = odeset('RelTol',10^(-tol),'AbsTol',10^(-tol));
    matlabRK = ode15s(f,tInt,y0,options);
    matlabWork15s(tol) = matlabRK.stats.nfevals+m*matlabRK.stats.npds;
    matlabError15s(tol) = norm(errorEstimator-matlabRK.y(:,end));
    matlabRK = ode23s(f,tInt,y0,options);
    matlabWork23s(tol) = matlabRK.stats.nfevals+m*matlabRK.stats.npds;
    matlabError23s(tol) = norm(errorEstimator-matlabRK.y(:,end));
end
fig = figure; %work vs error
loglog(ourWork,ourError,'*-b')
hold on
loglog(matlabWork15s,matlabError15s,'*-r')
loglog(matlabWork23s,matlabError23s,'*-g')
xlabel('Work = nfun + m \cdot njac')
ylabel('Error')
legend('RKs.m and onestep.m','ode15s','ode23s');
title('Work-precision diagram')
set(gca,'fontsize',15)
saveTightFigure(fig,'Figures/task10.pdf')
