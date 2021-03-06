function [tnext, ynext, le, iflag, nfun, njac] = onestep_Y3(f,jac,tn,yn,h,Tolit)
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

[A, c, g] = method();
m = length(yn);
maxIterations = 50;
%% To avoid errors:
tnext = 0;
ynext = 0;
le = 0;
iflag = -1;
nfun = 0;
njac = 1;%only one jacobian evaluation in this function
%% Newton iteration for finding stage values Y1, Y2, Y3, Y4
% Initializing Y to length of Y_0 vector
Y = zeros(m,4);
% Y_1 = y_n
Y(:,1) = yn;
%jac = double(jac(yn(1),yn(2),tn));
thisJac = jac(tn,yn);
I = eye(size(thisJac));
K = zeros(m,4);
K(:,1) = Y(:,1);
J = (I - h*g*thisJac);
%% Calculating Yi
for i = 2:3
    DY = Inf;
    iteration = 0;
% Finding Ki
    Y(:,i) = Y(:,i-1); %to minimize iterations in while loop
    K(:,i) = Y(:,1);
    for j = 1:(i-1)
        K(:,i) = K(:,i) + h*A(i,j)*f(tn + c(j)*h, Y(:,j));
        nfun = nfun+1;
    end
% Finding Yi nummericaly

        while double(norm(DY)) >= Tolit && iteration < maxIterations
            DY = J\(h*g*f(tn+c(i)*h,Y(:,i))-Y(:, i) + K(:,i));
            nfun = nfun+1;
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

tnext = tn + h;
ynext = Y(:, 3);

end