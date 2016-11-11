function jac = jacobian_real(name)
syms y1
syms y2
syms y3
syms t
if strcmp(name,'Linear test problem')
%Jacobi linear test problem
f1(y1,y2,t) = -2*y1+y2+t;
f2(y1,y2,t) = y1-2*y2+t+3;
J = [diff(f1,y1),diff(f1,y2);...
     diff(f2,y1),diff(f2,y2)];
 jac = J;%(a1, a2, time)
elseif strcmp(name,'Van der Pol equation')
%Jacobi Van der Pol equation
mu = 5;
x1(y1,y2) = y2;
x2(y1,y2) = mu*(1-y1^2)-y1;
JvP = [diff(x1,y1),diff(x1,y2);...
       diff(x2,y1),diff(x2,y2)];
   jac = JvP;
elseif strcmp(name,'The Robertson reaction')
%Jacobi Robertson reaction
z1 = -0.04*y1+10^4*y2*y3;
z2 = 0.04*y1-10^4*y2*y3-3*10^7*y2^2;
z3 = 3*10^7*y2^2;
JR = [diff(z1,y1),diff(z1,y2), diff(z1,y3);...
      diff(z2,y1),diff(z2,y2), diff(z2,y3);...
      diff(z3,y1),diff(z3,y2), diff(z3,y3)];
  jac = JR;
else
    fprintf('Ingen navn stemmer, bruk enten: \n test\n van der pol\n robertson\n')
end
end
  
  