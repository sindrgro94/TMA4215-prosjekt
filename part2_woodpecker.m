%Part 2 - The woodpecker toy problem
clear all;
close all;
clc;
a = 0.025; %m
m1 = 0.0003; %kg
d1 = 0.18335;
Ok1 = 10; %degrees
I2 = 7*10^-7; %kgm^2
b = 0.015; %m
m2 = 0.0045; %kg
d2 = 0.04766;
Ok2 = 12; %degrees
c = 0.0056;%Nm
g = 9.81; %m/s^2
q = I2+m2*b^2*(1-m2/(m1+m2));
%O0 = [
%phase a,c, d:
f1 = @(t,O) [O(2); (-c*O(1)+m2*b*g); 0; 0];
f2 = @(t,O) [O(2); (-c*O(1))/q; O(4); g+m2*b*c*O(1)/(q*m1+m2)];
Jac1 = @(t,O) [0 1 0 0; -c/(I2+m2*b^2) 0 0 0; 0 0 0 0; 0 0 0 0];
Jac2 = @(t,O) [0 1 0 0; -c/q 0 0 0; 0 0 0 1; m2*b*c/(q*(m1+m2)) 0 0 0];
