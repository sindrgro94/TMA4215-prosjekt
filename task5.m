%% Answering task 5
% Plotting stability domain [-a, a, -b, b]
% Heavily influenced by problem set 4
%% Arranging constant
clear all;
%syms z
a = 9; b = 9;
%a1 = 9; b1 = 6;
%a2 = 3; b2 = 3;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
%[x1, y1] = meshgrid(linspace(-a1, a1), linspace(-b1, b1));
%[x2, y2] = meshgrid(linspace(-a2, a2), linspace(-b2, b2));
z = x + i*y;
%z1 = x1 + i*y1;
%z2 = x2 + i*y2;

%% Stability function solving
%[gamma, R1Imported, R1HatImported] = task_2_to_4();
%R1FromFunction = abs(R1Imported);
%R1Substituted = subs(R1Imported, z, zDefinition);
%R1Double = double(R1Substituted);
%R1Real = real(R1Substituted);
%R1Plot = abs(R1Substituted);
%R1HatReal = abs(subs(R1Hat, z, zDefinition));
%% Explicit stability functions

R1Messy = (- (42748480185501301681749991164533*z.^3)./...
    2923003274661805836407369665432566039311865085952 +...
    (231375910929103659658712603652903*z.^2)./...
    162259276829213363391578010288128 + (16623663410109777*z)./...
    9007199254740992 - 6)./(6*((7851873215395081*z)./...
    18014398509481984 - 1).^3);
R1 = abs(R1Messy);

R1HatMessy = 1155326039345911*z./4503599627370496 - ...
    (58982356476257269084321494156719*z.^2)./...
    (162259276829213363391578010288128 + 2)./...
    (2*((7851873215395081*2)./18014398509481984 - 1).^2);
R1Hat = abs(R1HatMessy)
%% Making plot
%contourf(x, y, R1Real, [1 1], 'k')
%figure
%subplot(1,2,1)
fontSize = 12
%contourf(x1, y1, R1, [1 1], 'k')
contourf(x, y, R1, [1 1], 'k')
%axis equal, axis([-a1 a1 -b1 b1]), grid on
axis equal, axis([-a a -b b]), grid on
title('Plot of stability function $R(z)$', 'Interpreter', 'LaTeX',...
    'Fontsize', fontSize)
xlabel('$Re(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
ylabel('$Im(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
hold on
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
%plot([-a1, a1], [0, 0], 'k', 'LineWidth', 1);
%plot([0, 0], [-a1, a1], 'k', 'LineWidth', 1);
hold off

figure
%subplot(1,2,2)
contourf(x, y, R1Hat, [1 1], 'k')
%contourf(x2, y2, R1Hat, [1 1], 'k')
%axis equal, axis([-a2 a2 -b2 b2]), grid on
axis equal, axis([-a a -b b]), grid on
title('Plot of stability function $\widehat{R}(z)$', 'Interpreter', 'LaTeX',...
    'Fontsize', fontSize)
xlabel('$Re(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
ylabel('$Im(z)$','Interpreter', 'LaTeX','Fontsize', fontSize)
hold on
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
%plot([-a2, a2], [0, 0], 'k', 'LineWidth', 1);
%plot([0, 0], [-a2, a2], 'k', 'LineWidth', 1);
hold off



