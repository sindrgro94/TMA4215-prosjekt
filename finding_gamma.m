%finding gamma
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
yn_1 = subs(yn_1, l*h, z);
R(z) = yn_1/yn;
R_inf = limit(R(z),inf);


%ved å la z -> inf har jeg kommet frem til at
%R(inf) = -(6*g^3+18g^2-9g+1)/6g^3 hvilket gir
%g = 0.27


