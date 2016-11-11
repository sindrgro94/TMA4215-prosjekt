function onestep_solver(f,h,tn,yn,Tolit,testfunction)
%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @function jac jacobi matrix

%% Initializing constants and callings
clear all
N = ceil((tint(0)-tint(1))/h);
%jac = jacobi(parameter)

%% Initializing test functions
% All test problems are formatted {function, Jacobian string argument,
% current time step value, current solution value, stepsize and 

linearTestProblem = {f, ...
    testfunction, tint(0), yn};

testFunctions = {linearTestProblem};

%% Calling onestep function
% First establishing initial values
y1 = 1;
y2 = 2;
y3 = 0;
t = 0;
J = jacobian_real(testFunctions{1}{2}, y(1), y(2), y(3), tint(1));
tn = testFunctions{1}{3};
yn = testFunctions{1}{4};
iflag = 0;
y = zeros(2,N);
t = zeros(1,N);
% Looping trough N steps for finding the solution
tic
eg = 0; %Error globaly
for i = 1:N
    tic
    [tnext, ynext, le, iflag] = onestep(testFunctions{1}{1},...
        jac,tn,yn,h,...
        Tolit);
    toc
    try
        iflag = -1;
    catch
        warning('The Newton method ran out of iterations. Reduce stepsize and try again');
        tnext = 0;
        ynext = yn;
    end
    y(:,i) = ynext;
    t(i) = tnext;
    yn = ynext;
    tn = tnext;
    eg = eg +le;
end
toc
%% Plotting section






