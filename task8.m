function task8()
TestProblems = {'Linear test problem','Van der Pol equation','The Robertson reaction'};

h = 0.01;
tint = [0,1];
yn = [1;2];
f = {@(t,y) [t - 2*y(1) + y(2) ; t + y(1)- 2*y(2) + 3],...
     @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)],...
     @(t,y) [-0.04*y(1)+10^4*y(2)*y(3);...
             0.04*y(1)-10^4*y(2)*y(3)-3*10^7*y(2)^2;...
             3*10^7*y(2)^2]};

Tolit = 0.5;
onestep_solver(f{1},h,tint,yn,Tolit,TestProblems(1))




end