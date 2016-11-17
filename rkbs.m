function [t, y] = rkbs(f, t0, tend, y0, tol, h, P)
%
% Solve the ODE y' = f(t, y) with Bogacki-Shampines method. 
% Input: f       : Function handle to the rhs. of the ODE. 
%        t0, tend: The integration interval. 
%        y0      : Initial value.
%        tol     : Local error tolerance.
%        h       : Initial stepsize.
%        P       : Pessimist factor.
% Output: t      : Time steps.
%         y      : The solution at the points t.

% Rough preallocation of solution assuming constant stepsize.
N = (tend - t0) / h;
t = zeros(1, N);
y = zeros(length(y0), N);
t(1)    = t0;
y(:, 1) = y0;

% Current solution.
tn = t0;
yn = y0;

k1 = f(tn, yn); 

n = 2; % Not 1, because MATLAB starts indexing at 1.
while tn < tend
    % Make appropriate last step.
    if tn + h > tend
        h = tend - tn;
    end
    
    % One step.
    k2 = f(tn + .5  * h, yn + .5  * h * k1);
    k3 = f(tn + .75 * h, yn + .75 * h * k2);  
    y_new = yn + h * (2*k1 + 3*k2 + 4*k3) / 9;         % 3rd-order method.
    k4 = f(tn + h, y_new);
    z_new = yn + h * (7*k1 + 6*k2 + 8*k3 + 3*k4) / 24; % 2nd-order method.
    
    % Local error estimate.
    locerr_est = norm(y_new - z_new);
    if locerr_est < tol % The step is successful.
        yn = y_new;
        tn = tn + h;
        k1 = k4; % FSAL property.
    
        % Store the solution.
        t(n)    = tn;
        y(:, n) = yn;
        n = n + 1;
    end
    
    % New stepsize (EPS).
    h = P * (tol / locerr_est)^(1/3) * h;
end

% Remove possible unused slots.
if n-1 < N
    t(n:end)    = [];
    y(:, n:end) = [];
end
end