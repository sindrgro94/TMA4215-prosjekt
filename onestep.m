function [tnext, ynext, le, iflag] = onestep(f,jac,tn,yn,h,Tolit)
g = 0.435866762;
%for test example:
%f = @(t,y) [t - 2*y(1) + y(2);t + y(1)- 2*y(2) + 3];
%jac = jacobi('test')
%tn = 0
%yn = [1;2]
%h = 0.1
%Tolit = 0.1
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
%%
%g = sym('g'); %gamma
c = [0; 2*g; 1; 1];
bHat = [(-4*g^2+6*g-1)/(4*g); (-2*g+1)/(4*g); g; 0];
b = [(6*g-1)/(12*g); -1/(12*g*(2*g-1)); (-6*g^2+6*g-1)/(3*(2*g-1)); g];
A = [0,0,0,0;...
    g,g,0,0;...
    (-4*g^2+6*g-1)/(4*g), (-2*g+1)/(4*g), g, 0;...
    (6*g-1)/(12*g), -1/(12*g*(2*g-1)), (-6*g^2+6*g-1)/(3*(2*g-1)), g];
%% Newton iteration
Y = zeros(4,4)
Y(:,1) = yn';
%Y2
I = eye(size(jac));

%Y(:,2) = Y(:,1);
summen = 0
for i = 1:100
    for j = 1:2
        summen = summen + A(2,j)*f(tn+c(j)*h,Y(:,j));
    end
    DY = (I-h*g*jac)^-1*(-Y(2)+yn+h*summen);
    Y(:,2) = Y(:,2)+DY
end


















end