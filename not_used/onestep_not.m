function [tnext, ynext, le, iflag] = onestep(f,jac,tn,yn,h,Tolit)
% [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit)
% Do one step with an implicit RK?method method.

% Input arguments:
% f, jac: the functions f(t, y) and Jac(t, y);
% tn, yn: time and state variables
% h: step size
% Tolit: tolerance for the Newton iterations.

% Output arguments:
% tnext, ynext: time and state variables after one step
% le: Local error estimator.
% iflag = 1: Iterations are successfull
% = ?1: The iterations fails. t and y are not updated
%% Setting constants
%g = sym('g'); %gamma
g = 0.435866762;
c = [0; 2*g; 1; 1];
bHat = [(-4*g^2+6*g-1)/(4*g); (-2*g+1)/(4*g); g; 0];
b = [(6*g-1)/(12*g); -1/(12*g*(2*g-1)); (-6*g^2+6*g-1)/(3*(2*g-1)); g];
A = [0,0,0,0;...
    g,g,0,0;...
    (-4*g^2+6*g-1)/(4*g), (-2*g+1)/(4*g), g, 0;...
    (6*g-1)/(12*g), -1/(12*g*(2*g-1)), (-6*g^2+6*g-1)/(3*(2*g-1)), g];
m = length(yn);
maxIterations = 2;
%% Newton iteration for finding stage values Y1, Y2, Y3, Y4
% Initializing Y to length of Y_0 vector
Y = zeros(m,4);
% Y_1 = y_n
Y(:,1) = yn';
%Y2
I = eye(size(jac));
%Y(:,2) = Y(:,1);
K = zeros(m,4);
iteration = 0
DY = Inf
for i = 1:4
    for j = 1:(i - 1)
        K(:,i) = Y(:,1) + h*(A(i,j)*f(tn + c(j)*h, Y(:,j)))
    end
    while double(norm(DY)) >= Tolit && iteration < maxIterations
        %DY = inv(I - h*g*jac)*(-Y(:, j) - Y(:,i) + K(:, i))
        DY = inv(I - h*g*jac)*(- Y(:,i) + K(:, i))
        Y(:,i) = Y(:, i) + DY
        iteration = iteration + 1;
        if iteration == maxIterations
            % Testing whether k has ran out of maxIterations
            % If true, sets flag negative and returns to onestep function
            iflag = -1;
            tnext = NaN
            ynext = NaN
            le = NaN
            return
        end
        
    end
    if double(norm(DY)) < Tolit && iteration < maxIterations
        % Testing whether DY is less than given tolerance
        % If true, set returned iflag to 1 and break out of loop
        iflag = 1
        break
    end
    if i == maxIterations
        % Testing whether k has ran out of maxIterations
        % If true, sets flag negative and returns to onestep function
        iflag = -1;
        tnext = NaN
        ynext = NaN
        le = NaN
        return
    end
end

le = Y(:, 4) - Y(:, 3)
tnext = tn + h
ynext = Y(:, 4)


















end