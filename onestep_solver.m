%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @function jac jacobi matrix

%% Initializing constants and callings
clear all
N = [1000];
interval = 0;
%jac = jacobi(parameter)
t0 = 0;
tend = 1;
h = (tend - t0) / N;
Tolit = 0.5;

%% Initializing test functions
% All test problems are formatted {function, Jacobian string argument,
% current time step value, current solution value, stepsize and 

linearTestProblem = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3], ...
    'test', 0, [1 ; 2]};

testFunctions = {linearTestProblem};

%% Calling onestep function
% First establishing initial values
y1 = 1;
y2 = 2;
y3 = 0;
t = 0;
jac = jacobian_real(testFunctions{1}{2}, y1, y2, 0, t);
tn = testFunctions{1}{3};
yn = testFunctions{1}{4};
iflag = 0;
y = zeros(2,N);
t = zeros(1,N);
% Looping trough N steps for finding the solution
for i = 1:length(N)
    [tnext, ynext, le, iflag] = onestep(testFunctions{1}{1},...
        jac,tn,yn,h,...
        Tolit);
    try
        iflag = -1;
    catch
        warning('The Newton method ran out of iterations. Reduce stepsize and try again');
        tnext = 0;
        ynext = yn;
    end
    y(i) = ynext;
    t(i) = tnext;
    yn = ynext;
    tn = tnext;
end

%% Plotting section