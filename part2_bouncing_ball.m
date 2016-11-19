%Part 2 - Bouncing ball
clear all;
close all;
clc;
y0 = [10;0];
f = @(t,y) [y(2); -9.81];
%{use eventlocator, event, max stepsize, an event if y is ... than event}
eventLocator = {true,0,1e-3,'smaller'};
jac = @(t,y) [0 1; 0 0];
t0=0;
tend=100;
Tol=10^-8;
h0=10^-4;
for bounces = 1:40
    [t, y, iflag, nfun, njac] = RKs(f, jac, t0, tend, y0, Tol, h0,eventLocator);
    stop = find(y(1,:)<0);
    %find accurate t and y on impact:
    [tEvent,yEvent] = hermiteInterpolation(t((stop-1):stop),y(:,(stop-1):stop));
    %y0 = [0;-0.9*y(2,stop-1)];
    y0 = [0;-0.9*yEvent(2)];
    if bounces == 1
        ballY = [y(1,(1:stop-1)), yEvent(1)];
        ballT = [t(1:stop-1), tEvent];
    else
        ballY = [ballY, y(1,(1:stop-1)), yEvent(1)];
        ballT = [ballT, (ballT(end)+t(1:stop-1)), ballT(end)+tEvent];
    end 
end
h = figure;
hold on
plot(ballT,ballY)
xlabel('Time [s]')
ylabel('Height [m]')
set(gca,'fontsize',15)
hold off
saveTightFigure(h,'Figures/bouncing_ball.pdf')