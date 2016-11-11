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

[A, c, g, s, bHat, b] = method();
m = length(yn);
maxIterations = 3;

%% Newton iteration for finding stage values Y1, Y2, Y3, Y4
% Initializing Y to length of Y_0 vector
Y = zeros(m,4);
% Y_1 = y_n
Y(:,1) = yn';
I = eye(size(jac));
K = zeros(m,4);
K(:,1) = Y(:,1);
Ktest = zeros(m,4);
Ktest(:,1) = Y(:,1);
J = (I - h*g*jac);
%% Calculating Yi
for i = 2:4
    DY = Inf;
    iteration = 0;
% Finding Ki

    K(:,i) = Y(:,1);
    Ktest(:,i) = sum(h*A(i,1:i-1)*f(tn+c(1:i-1)*h,Y(:,1:i-1)))
    for j = 1:i-1
        K(:,i) = K(:,i) + h*A(i,j)*f(tn + c(j)*h, Y(:,j));
    end
    disp(K)
% Finding Yi nummericaly

        while double(norm(DY)) >= Tolit && iteration < maxIterations
            DY = J\(h*g*f(tn+c(i)*h,Y(:,i))-Y(:, i) + K(:,i));
            Y(:,i) = Y(:, i) + DY;
            iteration = iteration + 1;
        end

    if double(norm(DY)) < Tolit && iteration < maxIterations
        % Testing whether DY is less than given tolerance
        % If true, set returned iflag to 1 and break out of loop
        iflag = 1;
    end
    if iteration == maxIterations
        % Testing whether k has ran out of maxIterations
        % If true, sets flag negative and returns to onestep function
        iflag = -1;
        return
    end
end

le = abs(Y(:, 4) - Y(:, 3));
tnext = tn + h;
ynext = Y(:, 4);

end