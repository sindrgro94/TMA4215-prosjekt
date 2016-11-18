%Part 2 - The woodpecker toy problem
clear all;
close all;
clc;
a = 0.025; %m
m1 = 0.0003; %kg
d1 = 0.18335;
Ok1 = degtorad(10); %degrees
I2 = 7*10^-7; %kgm^2
b = 0.015; %m
m2 = 0.0045; %kg
d2 = 0.04766;
Ok2 = degtorad(12); %degrees
c = 0.0056;%Nm
g = 9.81; %m/s^2
q = I2+m2*b^2*(1-m2/(m1+m2));
maxStepSize = 10^-1;
%initial conditions:
t0 = 0;
tend = 10; %must be large enough.
Tol=10^-8;
h0=maxStepSize;
O0 = [degtorad(35);0;0;0];
%phase a,c, d:
f1 = @(t,O) [O(2); (-c*O(1)+m2*b*g)/(I2+m2*b^2); 0; 0];
Jac1 = @(t,O) [0 1 0 0; -c/(I2+m2*b^2) 0 0 0; 0 0 0 0; 0 0 0 0];
%phase b,e:
f2 = @(t,O) [O(2); (-c*O(1))/q; O(4); g+m2*b*c*O(1)/(q*m1+m2)];
Jac2 = @(t,O) [0 1 0 0; -c/q 0 0 0; 0 0 0 1; m2*b*c/(q*(m1+m2)) 0 0 0];
%{use eventlocator, event, max stepsize, an event if y is ... than event}
eventLocatorA = {true,Ok1,maxStepSize,'smaller'};
eventLocatorB = {true,-Ok1,maxStepSize,'smaller'};
eventLocatorC = {true,-Ok2,maxStepSize,'smaller'};
eventLocatorD = {true,-Ok1,maxStepSize,'bigger'};
eventLocatorE = {true,Ok1,maxStepSize,'bigger'};

for bounces = 1:5
    %%%%%%%%%%%STATE A:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorA);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorA{2}))
        return
    end
    stop = find(O(1,:)<=eventLocatorA{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorA{2});
    O0 = [eventLocatorA{2};OEvent(2:4)];
    %Update answer:
    if bounces == 1
        woodpeckerO = [O(:,(1:stop-1)), OEvent];
        woodpeckerT = [t(1:stop-1), tEvent];
    else
        woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
        woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
    end 
    
     %%%%%%%%%%%STATE B:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f2, Jac2, t0, tend, O0, Tol, h0,eventLocatorB);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorB{2}))
        return
    end
    stop = find(O(1,:)<=eventLocatorB{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorB{2});
    O0 = [eventLocatorB{2};OEvent(2:4)];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
    
     %%%%%%%%%%%STATE C:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorC);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorC{2}))
        return
    end
    stop = find(O(1,:)<=eventLocatorC{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorC{2});
    O0 = [eventLocatorC{2};OEvent(2:4)];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
    
     %%%%%%%%%%%STATE D:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorD);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorD{2}))
        return
    end
    stop = find(O(1,:)>=eventLocatorD{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorD{2});
    O0 = [eventLocatorD{2};OEvent(2:4)];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
    
     %%%%%%%%%%%STATE E:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f2, Jac2, t0, tend, O0, Tol, h0,eventLocatorE);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorE{2}))
        return
    end
    stop = find(O(1,:)>=eventLocatorE{2});
    %find accurate t and O on impact:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorE{2});
    O0 = [eventLocatorE{2};OEvent(2:4)];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
end