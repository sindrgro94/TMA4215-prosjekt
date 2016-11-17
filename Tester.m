%% Set up Lotka-Volterra system and plot solution and stepsizes.
f = @(t, y) [y(1) - y(1)*y(2); y(1)*y(2) - 2*y(2)];
y0 = [1; 2];
t0 = 0;
tend = 20;

[t, y] = rkbs(f, t0, tend, y0, 1e-5, 1/20, .9);
figure(1);
subplot(2, 1, 1)
plot(t, y)
xlim([t(1)-.5 t(end)+.5])
xlabel('$t$', 'Interpreter', 'LaTeX', 'Fontsize', 12);
l = legend('$y_1$', '$y_2$');
set(  l, 'Color', 'none', 'Location', 'northwest', 'Interpreter', 'LaTeX', 'Fontsize', 10);
set(gca, 'Color', 'none');

subplot(2, 1, 2)
plot(t(2:end), diff(t))
xlim([t(1)-.5 t(end)+.5])
xlabel(  '$t$', 'Interpreter', 'LaTeX', 'Fontsize', 12);
ylabel('$t_n$', 'Interpreter', 'LaTeX', 'Fontsize', 12);
set(gca, 'Color', 'none');

%% Plot of global error divided by tolerance.
tol = 2.^-(1:20);
globalerr = zeros(size(tol));

for i = 1:length(tol)
    [~, y] = rkbs(f, t0, tend, y0, tol(i), 1/20, .9);
    
    % Reference solution.
    options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
    [~, yref] = ode45(f, [t0, tend], y0, options);
    
    % Global error at the endpoint.
    globalerr(i) = norm(y(:, end) - yref(end, :)');
end

figure(2)
plot(globalerr ./ tol)
title('$\displaystyle\frac{\| e_N \|_2}{\texttt{TOL}}$', 'Interpreter', 'LaTeX', 'Fontsize', 12);
set(gca, 'Color', 'none');