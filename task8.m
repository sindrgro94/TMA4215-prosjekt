function task8()
%% 
TestProblems = {'Linear test problem','Van der Pol equation','The Robertson reaction'};

tint = [0,1];
mu = 50;
f = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3],...
     @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)],...
     @(t,y) [-0.04*y(1)+10^4*y(2)*y(3);...
             0.04*y(1)-10^4*y(2)*y(3)-3*10^7*y(2)^2;...
             3*10^7*y(2)^2]};
g = @(t,y) [t+exp(-t);t+exp(-t)+1];

Tolit = 10^-2;

%% Linear Test Equation
yn = [1,2];
h = zeros(2,2);
%eg = zeros(2,5);
cnt = 1;
for i = 1:0.5:3.5
    tic
    [~,y4,t,y3] = onestep_solver(f{1},10^-i,tint,yn,Tolit,TestProblems(1),mu);
    h(1,cnt) = 10^-i;
    opts = odeset('AbsTol', 10^(-12), 'RelTol', 10^(-12));
    %[~, yref] = ode15s(f{1}, tint, yn, opts); %solution to compare with
    yref = g(t(end),6);
    
    eg(1,cnt) = norm(y4(:, end) - yref);
    eg(2,cnt) = norm(y3(:, end) - yref);
    

    cnt = cnt+1;
    toc
end

figure();
loglog(h(1,:),eg)
hold on
loglog(h(1,:),h(1,:).^3/1000)
title('Global error estimate for the linear test problem')
legend('Numeric error Y_4','Numeric error Y_3','O(h^3)')
hold off
%% Van der Pol
cnt = 1;
yn = [2;0];
for i = 2:0.1:4
    tic
    [~,y4,t,y3] = onestep_solver(f{2},10^-i,tint,yn,Tolit,TestProblems(2),mu);
    h(2,cnt) = 10^-i;
    [~, yref] = ode15s(f{2}, [0,t(end)], yn, opts); %solution to compare with
    eg(1,cnt) = norm(y4(:, end) - yref(end,:)');
    eg(2,cnt) = norm(y3(:, end) - yref(end,:)');
    cnt = cnt+1;
    toc
end
loglog(h(2,:),eg)
hold on
loglog(h(2,:),h(2,:).^2)
title(sprintf('Global error for %s',TestProblems{2}))
legend(sprintf('Numeric error with mu = %i',mu),'O(h^2)')

%% Robertson
cnt = 1;
tint = [0,40];
yn = [1;0;0];
for i = 4%:0.5:5
    tic
    [eg(3,cnt),y,t] = onestep_solver(f{3},10^-i,tint,yn,Tolit,TestProblems(3),mu);
    h(3,cnt) = 10^-i;
    cnt = cnt+1;
    toc
end
figure();
loglog(h(3,:),eg(3,:))
hold on
loglog(h(3,:),h(3,:).^2)
title('Global error for Robertson equation')
legend('Numeric error','O(h^2)')
    





end