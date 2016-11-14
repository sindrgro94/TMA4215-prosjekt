function jac = jacobi(name)
if strcmp(name,'Linear test problem')
%Jacobi linear test problem
jac = @(y,t,mu) [-2 1; 1 -2];
elseif strcmp(name,'Van der Pol equation')
%Jacobi Van der Pol equation
jac = @(y,t,mu) [0, 1; -10*y(1)-1 , 0];
elseif strcmp(name,'The Robertson reaction')
%Jacobi Robertson reaction
jac = @(y,t,mu) [-1/25, 10000*y(3), 10000*y(2);...
                 1/25, -60000000*y(2) - 10000*y(3), -10000*y(2);...
                 0, 60000000*y(2),0];
else
    fprintf('Ingen navn stemmer, bruk enten: \n Linear test problem\n')
    fprintf('Linear test problem\n The Robertson reaction\n')
end
end
  
  
  
  