%% Answering task 5
% Plotting stability domain [-a, a, -b, b]
% Heavily influenced by problem set 5
%% Arranging constant
clear all;
a = 4; b = 4;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
z = x + i*y;

%% Stability function solving
[gamma, R1, R1_hat] = task_2_to_4()
stabilityFunction = subs(R1, g, gamma)
stabilityFunctionHat = subs(R1_hat, g, gamma)



