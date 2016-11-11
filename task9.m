% Task 9
clear all;
clc;
TestProblems = {'Linear test problem','Van der Pol equation','The Robertson reaction'};
mu = 10; %can be chosen
f = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3],...
     @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)],...
     @(t,y) [-0.04*y(1)+10^4*y(2)*y(3);...
             0.04*y(1)-10^4*y(2)*y(3)-3*10^7*y(2)^2;...
             3*10^7*y(2)^2]};
jac = {0,0,0};
for i = 1:3
    jac{i} = jacobian_real(TestProblems{i});
end
tInt = {[0 1],[0 mu],[0 40]};
y0 = {[1;2],[2;0],[1,0,0]};
Tol = {10^(-2), 10^(-2), 10^(-2)}; %can be chosen
h0 = 0.1; %can be chosen
for i = 1:3
    [t, y, iflag, nfun, njac] = RKs(f{i}, jac{i}, tInt{i}(1), tInt{i}(2), y0{i}, Tol{i}, h0);
    if iflag == 1 %...we make plots:
        stop = find(t == tInt{i}(2));
        figure;
        title(TestProblems{i});
        hold on
        subplot(2,1,1)
        plot(t(1:stop),y(1:stop));
        subplot(2,1,2)
        plot(t(1:stop-1),(t(2:stop)-t(1:stop-1))); 
    else
        fprintf(TestProblems{i});
        fprintf(' went wrong.\n');
    end
end
