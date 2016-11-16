%Part 2 - Bouncing ball
clear all;
close all;
clc;
y0 = [10;0];
f = @(t,y) [y(2); -9.81];
eventLocator = {true,0,1e-3};
jac = @(t,y) [0 1; 0 0];
t0=0;
tend=100;
Tol=10^-8;
h0=10^-2;
for bounces = 1:20
    [t, y, iflag, nfun, njac] = RKs(f, jac, t0, tend, y0, Tol, h0,eventLocator);
    stop = find(y(1,:)<0);
    %find accurate t and y on impact:
    [tImpact,yImpact,speedImpact] = hermiteInterpolation(t((stop-2):stop),y(:,(stop-2):stop));
    %y0 = [0;-0.9*y(2,stop-1)];
    y0 = [0;-0.9*speedImpact];
    if bounces == 1
        ballY = [y(1,(1:stop-1)), yImpact];
        ballT = [t(1:stop-1), tImpact];
    else
        ballY = [ballY, y(1,(1:stop-1)), yImpact];
        ballT = [ballT, (ballT(end)+t(1:stop-1)), tImpact];
    end 
end
plot(ballT,ballY)