%% Matrise
%%%oppgave 1-2

%% setter konstanter
g = sym('g','real'); %gamma
syms h
c = [0; 2*g; 1; 1];
bHat = [(-4*g^2+6*g-1)/(4*g); (-2*g+1)/(4*g); g; 0];
b = [(6*g-1)/(12*g); -1/(12*g*(2*g-1)); (-6*g^2+6*g-1)/(3*(2*g-1)); g];
A = [0,0,0,0;...
    g,g,0,0;...
    (-4*g^2+6*g-1)/(4*g), (-2*g+1)/(4*g), g, 0;...
    (6*g-1)/(12*g), -1/(12*g*(2*g-1)), (-6*g^2+6*g-1)/(3*(2*g-1)), g];
I = eye(4);
one = ones(4);
%% Regner ut stabilitetsfunksjoner
R1intermediate = b' * inv(I - h*A)*one;
R2intermediate = bHat' * inv(I - h*A)*one;
R1 = simplify(sum(R1intermediate))
R2 = simplify(sum(R2intermediate))

%% Skrapseksjon
%term1 = inv(I - h*A)
%result = b' * term1
%jiha = 1 + term1
%enkel = simplify(sum(jiha))
