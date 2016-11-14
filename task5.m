%% Answering task 5
% Plotting stability domain [-a, a, -b, b]
% Heavily influenced by problem set 4
%% Arranging constant
clear all;
%syms z
a = 0.5; b = 0.5;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
z = x + i*y;

%% Stability function solving
%[gamma, R1, R1Hat] = task_2_to_4();
%R1Substituted = subs(R1, z, zDefinition);
%R1Double = double(R1Substituted);
%R1Real = real(R1Substituted);
%R1Plot = abs(R1Substituted);
%R1HatReal = abs(subs(R1Hat, z, zDefinition));
%% Explicit stability functions
R1Piggy = (- (42748480185501301681749991164533*z^3)/2923003274661805836407369665432566039311865085952 + (231375910929103659658712603652903*z^2)/162259276829213363391578010288128 + (16623663410109777*z)/9007199254740992 - 6)/(6*((7851873215395081*z)/18014398509481984 - 1)^3);
R1 = abs(R1Piggy)
%% Making plot
%contourf(x, y, R1Real, [1 1], 'k')
figure
contourf(x, y, R1, [1 1], 'k')
axis equal, axis([-a a -b b]), grid on
hold on
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
hold off



