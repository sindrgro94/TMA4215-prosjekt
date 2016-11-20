%% Answering task 5
% Plotting stability domain [-a, a, -b, b]
% Heavily influenced by problem set 4
%% Arranging constant
clear all;
%syms z
a = 9; b = 9;
%a1 = 9; b1 = 6;
a2 = 100; b2 = 100;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
%[x1, y1] = meshgrid(linspace(-a1, a1), linspace(-b1, b1));
[x2, y2] = meshgrid(linspace(-a2, a2), linspace(-b2, b2));
z = x + i*y;
%z1 = x1 + i*y1;
z2 = x2 + i*y2;
%% Explicit stability functions

R1Messy = (- (42748480185501301681749991164533*z.^3)./...
    2923003274661805836407369665432566039311865085952 +...
    (231375910929103659658712603652903*z.^2)./...
    162259276829213363391578010288128 + (16623663410109777*z)./...
    9007199254740992 - 6)./(6*((7851873215395081*z)./...
    18014398509481984 - 1).^3);
R1 = abs(R1Messy);

%R1HatMessy = 1155326039345911*z./4503599627370496 -(58982356476257269084321494156719*z.^2)./ (162259276829213363391578010288128 + 2)./ (2*((7851873215395081*2)./18014398509481984 - 1).^2);
R1HatMessy = -(58982356476257269084321494156719*z2.^2 - 41625007362317406238789597134848*z2 - 324518553658426726783156020576256)./(7851873215395081*z2 -18014398509481984).^2
R1Hat = abs(R1HatMessy) 
%% Plotting R(z)
%contourf(x, y, R1Real, [1 1], 'k')
fig = figure;
subplot(1,2,1)
fontSize = 15
%contourf(x1, y1, R1, [1 1], 'k')
contourf(x, y, R1, [1 1], 'k')
%axis equal, axis([-a1 a1 -b1 b1]), grid on
axis equal, axis([-a a -b b]), grid on
title('Stability function $R(z)$', 'Interpreter', 'LaTeX',...
    'Fontsize', fontSize)
xlabel('$Im(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
ylabel('$Re(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
hold on
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
%plot([-a1, a1], [0, 0], 'k', 'LineWidth', 1);
%plot([0, 0], [-a1, a1], 'k', 'LineWidth', 1);
set(gca,'fontsize',12)
hold off
% Plotting RHat(z)
%figure
h = subplot(1,2,2)
%contourf(x, y, R1Hat, [1 1], 'k')
contourf(x2, y2, R1Hat, [1 1], 'k')
axis equal, axis([-a2 a2 -b2 b2]), grid on
%axis equal, axis([-a a -b b]), grid on
title('Stability function $\widehat{R}(z)$', 'Interpreter', 'LaTeX',...
    'Fontsize', fontSize)
xlabel('$Im(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
ylabel('$Re(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
hold on
%plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
%plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
plot([-a2, a2], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a2, a2], 'k', 'LineWidth', 1);
set(gca,'fontsize',12)
%% Extremely neat method for adjusting subplot parameters
p = get(h, 'pos');
p(1) = p(1) + 0.1;
set(h, 'pos', p);
hold off
saveTightFigure(fig,'Figures/stabPlot.pdf')


