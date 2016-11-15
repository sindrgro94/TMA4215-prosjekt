function [eg,y,t] = onestep_solver(f,h,tint,yn,Tolit,testfunction,mu)
%% Introduction
% This script answers task 8
% For each iteration, calls onestep.m to solve one step of a test problem
%' @function jac jacobi matrix

%% Initializing constants and callings
N = ceil((tint(2)-tint(1))/h);
%% Finner Jacobi
jac = jacobi(testfunction{1},mu); %Her vil du f? riktig jacobi, men det er 
    %en funksjon av t og y, s? du m? gi inn t og hele y.
% if strcmp(testfunction{1},'Linear test problem')
%     jac = @(y) [-2,1;1,-2];
% elseif strcmp(testfunction{1},'Van der Pol equation')
%     jac = @(y) [0,1;-10*y(1)-1,0];
% elseif strcmp(testfunction{1},'The Robertson reaction')
%     jac = @(y) [-1/25,  10000*y(3),                     10000*y(2);...
%                 1/25,   - 60000000*y(2) - 10000*y(3),   -10000*y(2);...
%                 0,      60000000*y(2),                   0];
% else
%     fprintf('Ikke noe navnt registrert')
%     return
% end
%% Calling one step function
iflag = 0;
y = zeros(2,N);
y(:,1) = yn;
t = zeros(1,N);
tn = t(1);
% Looping trough N steps for finding the solution
eg = 0; %Error globaly
for i = 2:N
    %J = double(jac(yn(1),yn(2),tn)); % Slik at vi sender inn en matrise av doubles
    [tnext, ynext, le, iflag] = onestep(f,jac(yn(1)),t(i-1),y(:,i-1),h,Tolit);
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
%% Plotting section






