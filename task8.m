function task8()
%% Beskrivelse av problemer
% Jeg har jobbet litt med dette scriptet i kveld 13.11.16. Jeg er ganske
% sikker p? at vi kan effektivisere det en del ved ? endre m?ten vi finner
% jacobien p?. Er mulighens best om vi lagerer den p? en annen m?te slik
% Sindre snakket om, eller om vi faktisk kan approksimere den slik det st?r i boken
% side 286. Tror kanskje sistnevnte blir litt for un?yaktig... Jeg la til
% en linje (37) i onestep.m, grunnen er kommentert ved siden av.
% Totaltsett, med symbolsk jacobi bruker vi ca 0.01 sek per nye y-verdi, for symbolsk jacobi, som
% egt er ganske mye. Med en jacobi funksjon som i linje 9-21 i
% onestep_solver.m bruker den rundt 0.00007 sek per iterasjon. S? jeg
% stemmer egt for den jeg da.

% Problemene n? er at den ikke fungerer generelt for funksjonene v?re, men jeg tror kunne noen sm? justeringer
% kan fikse det

%% 
TestProblems = {'Linear test problem','Van der Pol equation','The Robertson reaction'};

%h = 0.01;
tint = [0,1];
mu = 50;
f = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3],...
     @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)],...
     @(t,y) [-0.04*y(1)+10^4*y(2)*y(3);...
             0.04*y(1)-10^4*y(2)*y(3)-3*10^7*y(2)^2;...
             3*10^7*y(2)^2]};

Tolit = 10^-3;

%% Linear Test Equation
yn = [1,2];
h = zeros(2,5);
eg = zeros(2,5);
cnt = 1;
for i = 2:0.5:4
    tic
    [eg(1,cnt),y,t] = onestep_solver(f{1},10^-i,tint,yn,Tolit,TestProblems(1),mu);
    h(1,cnt) = 10^-i;
    cnt = cnt+1;
    toc
end
figure();
loglog(h(1,:),eg(1,:))
hold on
loglog(h(1,:),h(1,:).^2)
title('Global error estimate for the linear test problem')
legend('Numeric error','O(h^2)')
%% Van der Pol
cnt = 1;
yn = [2;0];
for i = 3:0.5:5
    tic
    [eg(2,cnt),y,t] = onestep_solver(f{2},10^-i,tint,yn,Tolit,TestProblems(2),mu);
    h(2,cnt) = 10^-i;
    cnt = cnt+1;
    toc
end
figure();
loglog(h(2,:),eg(2,:))
hold on
loglog(h(2,:),h(2,:).^2)
title(sprintf('Global error for %s',TestProblems{2}))
legend(sprintf('Numeric error with mu = %i',mu),'O(h^2)')

%% Robertson
cnt = 1;
tint = [0,40];
yn = [1;0;0];
for i = 4:5
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