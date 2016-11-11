function onestep_solver(f,h,tint,yn,Tolit,testfunction)
%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @function jac jacobi matrix

%% Initializing constants and callings
N = ceil((tint(2)-tint(1))/h);
%% Calling onestep function
% First establishing initial values
t = 0;
jac = jacobian_real(testfunction);
iflag = 0;
y = zeros(2,N);
y(:,1) = yn;
t = zeros(1,N);
% Looping trough N steps for finding the solution
tic
eg = 0; %Error globaly
for i = 2:N
    tic
    [tnext, ynext, le, iflag] = onestep(f,jac,t(i-1),y(:,i-1),h,Tolit);
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






