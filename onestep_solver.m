%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @jac jacobi matrix

%% initializing test functions
% All test problems are formatted {function, Jacobian, current time, .....
% current solution value, stepsize and 

jac = jacobi('test')

f = @(t,y) [t - 2*y(1) + y(2);t + y(1)- 2*y(2) + 3];

[tnext, ynext, le, iflag] = onestep(f,jac,tn,yn,h,Tolit)