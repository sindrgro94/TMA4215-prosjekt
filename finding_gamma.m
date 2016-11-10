function x = finding_gamma()

%% Establishes the butcher tableu 
g = sym('g','real'); %gamma
c = [0; 2*g; 1; 1];
bHat = [(-4*g^2+6*g-1)/(4*g); (-2*g+1)/(4*g); g; 0];
b = [(6*g-1)/(12*g); -1/(12*g*(2*g-1)); (-6*g^2+6*g-1)/(3*(2*g-1)); g];
A = [0,0,0,0;...
    g,g,0,0;...
    (-4*g^2+6*g-1)/(4*g), (-2*g+1)/(4*g), g, 0;...
    (6*g-1)/(12*g), -1/(12*g*(2*g-1)), (-6*g^2+6*g-1)/(3*(2*g-1)), g];
%% RK4
% A = [0, 0, 0, 0;...
%     1/2,0, 0, 0;...
%     0, 1/2,0, 0;...
%     0,  0, 1, 0];
% b = [1/6,1/3,1/3,1/6];
% c = [0,1/2,1/2,1];
%% finding ki for test problem
yn = sym('yn');
l = sym('l');
h = sym('h');
z = sym('z');
%k1 = sym('k1'); k2 = sym('k2'); k3 = sym('k3'); k4 = sym('k4');
k1 = l*yn;
%k2 = l*(yn+h*g*k1+h*g*k2); løser denne for k2
k2 = l*yn*(1+h*A(2,1)*l)/(1-h*A(2,2)*l);
k3 = l*(yn+h*(A(3,1)*k1+A(3,2)*k2))/(1-h*l*A(3,3));
k4 = l*(yn+h*(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3))/(1-h*l*A(4,4));
%% Creating stability function R(z)
yn_1= simplify(yn+h*(b(1)*k1+b(2)*k2+b(3)*k3+b(4)*k4));
yn_1h = simplify(yn+h*(bHat(1)*k1+bHat(2)*k2+bHat(3)*k3+bHat(4)*k4));
yn_1 = subs(yn_1, l*h, z);
yn_1h = subs(yn_1h,l*h,z);
R(z,g) = yn_1/yn;
Rh(z,g) = yn_1h/yn;
R_inf(g) = limit(R(z,g),inf);
R_infh(g) = limit(Rh(z,g),inf);
%% Newtons method
x = 0.4;
differ = 0.4;
dR_inf = diff(R_inf,g);
while abs(differ) > (10^4)*eps
    x_np = x-R_inf(x)/dR_inf(x);
    differ = x_np-x;
    x = x_np;
end
x = double(x);
%% R_hat solution
y = solve(R_infh,g);
y_ret = y(1);
end

