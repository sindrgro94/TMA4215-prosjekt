% TMA4215 prosjekt 
% Oppgave 1-2

%% setter konstanter
%g = sym('g','real');    % gamma
syms h g  z               % gjer det samme som sym('variabel')?
c = [0; 2*g; 1; 1];
%% 
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
fprintf('\n\n')
fprintf('Utregna og forenkla stabilitetsfunksjon R1:\n\n')
R1 = simplify(sum(R1intermediate))
fprintf('Utregna og forenkla stabilitetsfunksjon R1:\n\n')
R2 = simplify(sum(R2intermediate))
R1z = subs(R1, g*h, z)
R2z = subs(R2, g*h, z)

%% Regner ut stabilitetskonstanter
R1Inf = limit(R1z, z, Inf)
R2Inf = limit(R2z, z, Inf)

%% Om A var invertibel /ikke-singulær
%R1InfI = 1 - b'*inv(A)*one
[V,D] = eig(A)

%% Newtons metode for å regne ut gamma
g = 0.05;
g_old = 100;
%g_true = 0.0623776;
iter = 0;
while abs(g_old-g) > 10^-3 && g ~= 0
    g_old = g;
    g = g - (g^3 - 0.165*g^2 + 3.993*10^-4)/(3*g^2 - 0.33*g);
    iter = iter + 1;
    fprintf('Iteration %d: x=%.20f, err=%.20f\n', iter, g, g_true-g);
    pause;
end

%% Skrapseksjon
%term1 = inv(I - h*A)
%result = b' * term1
%jiha = 1 + term1
%enkel = simplify(sum(jiha))
