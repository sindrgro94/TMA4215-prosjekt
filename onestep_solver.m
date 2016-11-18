function [y4,y3,t] = onestep_solver(f,h,tint,yn,Tolit,testfunction,mu)
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
y4 = zeros(length(yn),N);
y4(:,1) = yn;
y3 = zeros(length(yn),N);
y3(:,1) = yn;
t = linspace(tint(1),tint(2),N);
for i = 2:N
    [t(i), y4(:,i), ~, iflag,~,~] = onestep(f,jac,t(i-1),y4(:,i-1),h,Tolit);
    if iflag == -1;
        warning('The Newton method for onestep.m ran out of iterations. Reduce stepsize and try again');
        t(i) = 0;
        return
    end
    [~,y3(:,i),~,iflag,~,~] = onestep_Y3(f,jac,t(i-1),y3(:,i-1),h,Tolit);
    if iflag == -1;
        warning('The Newton method for onestep_Y3.m ran out of iterations. Reduce stepsize and try again');
        t(i) = 0;
        return
    end
    
end
end

