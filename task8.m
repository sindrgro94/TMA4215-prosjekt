function task8()
%% Initiating
TestProblems = {'Linear test problem','Van der Pol equation','The Robertson reaction'};

tint = [0,1];
mu = 50;
f = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3],...
     @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)],...
     @(t,y) [-0.04*y(1)+10^4*y(2)*y(3);...
             0.04*y(1)-10^4*y(2)*y(3)-3*10^7*y(2)^2;...
             3*10^7*y(2)^2]};
g = @(t,y) [t+exp(-t);t+exp(-t)+1]; %Analytic solution for the test problem

Tolit = 10^-8;

%% Linear Test Equation
yn = [1,2];
h = zeros(1,2);
cnt = 1;
for i = 0.5:0.5:3.5
    tic
    h(cnt) = 10^-i;
    [y4,y3,t] = onestep_solver(f{1},h(cnt),tint,yn,Tolit,TestProblems(1),mu);
    opts = odeset('AbsTol', 10^(-12), 'RelTol', 10^(-12));
    yref = g(t(end),6);
    eg(1,cnt) = norm(y4(:, end) - yref);
    eg(2,cnt) = norm(y3(:, end) - yref);
    cnt = cnt+1;
    toc
end

figure();
loglog(h(1,:),eg)
hold on
loglog(h, h.^2/100)
loglog(h, h.^3/1000)
title('Global error estimate for the linear test problem')
legend('Global error for advancing method','Global error for error estimating method','O(h^2)','O(h^3)')
hold off

ordY4 = polyfit(log(h),log(eg(1,:)),1);
ordY3 = polyfit(log(h),log(eg(2,:)),1);

orders(1,:) = [ordY4(1),ordY3(1)];
%% Van der Pol
cnt = 1;
yn = [2;0];
eg = zeros(2,4);
h = [0];
for i = 2:0.1:3.5
    tic
    [y4,y3,t] = onestep_solver(f{2},10^-i,tint,yn,Tolit,TestProblems(2),mu);
    h(cnt) = 10^-i;
    [~, yref] = ode15s(f{2}, [0,t(end)], yn, opts); %solution to compare with
    eg(1,cnt) = norm(y4(:, end) - yref(end,:)');
    eg(2,cnt) = norm(y3(:, end) - yref(end,:)');
    cnt = cnt+1;
    toc
end
figure();
loglog(h,eg)
hold on
loglog(h,h.^2/50000)
loglog(h,h.^3/1000)
title(sprintf('Global error for %s',TestProblems{2}))
legend('Global error for advancing method','Global error for error estimating method','O(h^2)','O(h^3)')
hold off

ordY4 = polyfit(log(h),log(eg(1,:)),1);
ordY3 = polyfit(log(h),log(eg(2,:)),1);

orders(2,:) = [ordY4(1),ordY3(1)];

%% Robertson
cnt = 1;
yn = [1;0;0];
eg = [0;0];
h = [0];
for i = 3:0.1:4
    tic
    h(cnt) = 10^-i;
    [y4,y3,t] = onestep_solver(f{3},h(cnt),tint,yn,Tolit,TestProblems(3),mu);
    [~,yref] = ode15s(f{3}, [0,t(end)], yn, opts);
    eg(1,cnt) = norm(y4(:, end) - yref(end,:)');
    eg(2,cnt) = norm(y3(:, end) - yref(end,:)');
    cnt = cnt+1;
    toc
end
figure();
loglog(h,eg)
hold on
loglog(h,h.^2/1000)
loglog(h,h.^3/10)
title(sprintf('Global error for %s',TestProblems{2}))
legend('Global error for advancing method','Global error for error estimating method','O(h^2)','O(h^3)')
hold off

ordY4 = polyfit(log(h),log(eg(1,:)),1);
ordY3 = polyfit(log(h),log(eg(2,:)),1);

orders(3,:) = [ordY4(1),ordY3(1)];

%% Conclusions
%Numerical orders displayed in the command window
disp(orders)

end