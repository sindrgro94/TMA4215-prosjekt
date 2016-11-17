function [y4,y3,t,eg] = onestep_solver(f,h,tint,yn,Tolit,testfunction,mu)
%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @function jac jacobi matrix

%% Initializing constants and callings
N = ceil((tint(2)-tint(1))/h)+1;
%% Finner Jacobi
jac = jacobi(testfunction{1},mu); %Her vil du f? riktig jacobi, men det er 
    %en funksjon av t og y, s? du m? gi inn t og hele y.

%% Calling one step function
iflag = 0;
y4 = zeros(length(yn),N);
y4(:,1) = yn;
y3 = zeros(length(yn),N);
y3(:,1) = yn;
t = linspace(tint(1),tint(2),N);
eg = 0;
for i = 2:N
    [t(i), y4(:,i), le, iflag,~,~,y3(:,i)] = onestep(f,jac,t(i-1),y4(:,i-1),h,Tolit);
    try
        iflag = -1;
    catch
        warning('The Newton method ran out of iterations. Reduce stepsize and try again');
        t(i) = 0;
    end
    eg = eg+le;
end
%% Plotting section




end

