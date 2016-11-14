% Task 10 - work-precision Diagram: (for Van der Pol equation)
mu = 10;
f = @(t,y) [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
jac = jacobi('Van der Pol equation');
tInt = [0 mu];
y0 = [2;0];
m = 2; %from y
h0 = 0.1; %can be chosen
for Tol = 2:10
    [t, y, iflag, nfun, njac] = RKs(f, jac, tInt(1), tInt(2), y0, 10^(-Tol), h0, mu);
    work = nfun+m*njac;
    sol = ode23s(f,tInt,y0,options);
    
end
%ode15s til eksakt l?sning 10^-12 feil