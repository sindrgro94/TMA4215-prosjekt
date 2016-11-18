function jac = jacobi(name,mu)
if strcmp(name,'Linear test problem')
%Jacobi linear test problem
jac = @(t,y) [-2 1; 1 -2];
elseif strcmp(name,'Van der Pol equation')
%Jacobi Van der Pol equation
jac = @(t,y) [0, 1; -2*mu*y(1)*y(2)-1 , mu*(1-y(1)^2)];
elseif strcmp(name,'The Robertson reaction')
%Jacobi Robertson reaction
jac = @(t,y) [-1/25, 10000*y(3), 10000*y(2);...
                 1/25, -60000000*y(2) - 10000*y(3), -10000*y(2);...
                 0, 60000000*y(2),0];
else
    fprintf('Ingen navn stemmer, bruk enten: \n Linear test problem\n')
    fprintf('Linear test problem\n The Robertson reaction\n')
end
end
  
  
  
  