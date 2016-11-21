%Part 2 - The woodpecker toy problem
clear all;
close all;
clc;
[a,m1,d1,Ok1,I2,b,m2,d2,Ok2,c,g] = constants();
q = I2+m2*b^2*(1-m2/(m1+m2));
g = -g;
maxStepSize = 10^-1;
%initial conditions:
t0 = 0;
tend = 1; %must be large enough.
Tol=10^-8;
h0=maxStepSize;
O0 = [degtorad(10);25;0;0];
%phase a,c, d:
f1 = @(t,O) [O(2); (-c*O(1)+m2*b*g)/(I2+m2*b^2); 0; 0];
Jac1 = @(t,O) [0 1 0 0; -c/(I2+m2*b^2) 0 0 0; 0 0 0 0; 0 0 0 0];
%phase b,e:
f2 = @(t,O) [O(2); (-c*O(1))/q; O(4); g + m2*b*c*O(1)/(q*(m1+m2))];
Jac2 = @(t,O) [0 1 0 0; -c/q 0 0 0; 0 0 0 1; m2*b*c/(q*(m1+m2)) 0 0 0];
%{use eventlocator, event, max stepsize, an event if y is ... than event}
eventLocatorA = {true,Ok1,maxStepSize,'smaller'};
eventLocatorB = {true,-Ok1,maxStepSize,'smaller'};
eventLocatorC = {true,-Ok2,maxStepSize,'smaller'};
eventLocatorD = {true,-Ok1,maxStepSize,'bigger'};
eventLocatorE = {true,Ok1,maxStepSize,'bigger'};
%impact sleve and beak:
impactSleveTop = @(theta_,z_) (1-d2)*(theta_+m2*b/(I2+m2*b^2)*(-z_));
impactSleveBottom = @(theta_,z_) (1-d1)*(theta_+((m2*b)/(I2+m2*b^2)*(-z_)));
impactBeak = @(theta_) -theta_;
%modified model:
% impactSleveTop = @(theta_,z_) (1-d1)*(theta_+m2*b/(I2+m2*b^2)*z_);
% impactSleveBottom = @(theta_,z_) (1-d2)*(theta_+((m2*b)/(I2+m2*b^2)*z_));

for bounces = 1:15
    disp(bounces)
    %%%%%%%%%%%STATE A:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorA);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorA{2}))
        return
    end
    stop = find(O(1,:)<eventLocatorA{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorA{2});
    %new initial conditions for state b:
    O0 = [eventLocatorA{2}; OEvent(2:4)];
    %Update answer:
    if bounces == 1
        woodpeckerO = [O(:,(1:stop-1)), OEvent];
        woodpeckerT = [t(1:stop-1), tEvent];
    else
        woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
        woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
    end 
%plot(radtodeg(O(1,(1:stop-1))),O(2,(1:stop-1)))
%hold on
     %%%%%%%%%%%STATE B:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f2, Jac2, t0, tend, O0, Tol, h0,eventLocatorB);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorB{2}))
        return
    end
    stop = find(O(1,:)<eventLocatorB{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorB{2});
    %new initial conditions for state c(z speed is 0):
    O0 = [eventLocatorB{2}; impactSleveTop(OEvent(2),OEvent(4)); OEvent(3); 0];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
     %%%%%%%%%%%STATE C:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorC);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorC{2}))
        return
    end
    if t(end) ~= tend %in case it did not reach the pole
        stop = find(O(1,:)<eventLocatorC{2});
    else
        woodpeckerO = [woodpeckerO, O(:,:)];
        woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t)];
        fprintf('The woodpecker did not reach the pole\n');
        return
    end
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorC{2});
    %new initial conditions for state d:
    O0 = [eventLocatorC{2};impactBeak(OEvent(2)); OEvent(3);0];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
     %%%%%%%%%%%STATE D:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f1, Jac1, t0, tend, O0, Tol, h0,eventLocatorD);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorD{2}))
        return
    end
    stop = find(O(1,:)>eventLocatorD{2});
    %find accurate t and O on event:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorD{2});
    %new initial conditions for state e:
    O0 = [eventLocatorD{2}; OEvent(2:3);0];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
     %%%%%%%%%%%STATE E:%%%%%%%%%%%%%%%%
    [t, O, iflag] = RKs(f2, Jac2, t0, tend, O0, Tol, h0,eventLocatorE);
    if iflag == -1
        fprintf('Error in RKs at theta K = %d\n',radtodeg(eventLocatorE{2}))
        return
    end
    stop = find(O(1,:)>eventLocatorE{2});
    %find accurate t and O on impact:
    [tEvent,OEvent] = hermiteInterpolationPecker(t((stop-1):stop),O(:,(stop-1):stop),eventLocatorE{2});
    %new initial conditions for state a(z speed is 0):
    O0 = [eventLocatorE{2};impactSleveBottom(OEvent(2),OEvent(4)); OEvent(3); 0];
    %Update answer:
    woodpeckerO = [woodpeckerO, O(:,(1:stop-1)), OEvent];
    woodpeckerT = [woodpeckerT, (woodpeckerT(end)+t(1:stop-1)), woodpeckerT(end)+tEvent];
end
angleVel = figure;
hold on
plot(radtodeg(woodpeckerO(1,:)),woodpeckerO(2,:))
xlabel('Angle [Degrees]')
ylabel('Angular Velocity [rad/s]')
set(gca,'fontsize',18)
hold off
heightTime = figure;
hold on
plot(woodpeckerT,woodpeckerO(3,:))
ylabel('Height [m]')
xlabel('Time [s]')
set(gca,'fontsize',18)
hold off
angleTime = figure;
hold on
plot(woodpeckerT,radtodeg(woodpeckerO(1,:)))
xlabel('Time [s]')
ylabel('Angle [Degrees]')
set(gca,'fontsize',18)
hold off
% saveTightFigure(angleVel,'Figures/angleVelocity_modified.pdf')
% saveTightFigure(heightTime,'Figures/heightTime_modified.pdf')
% saveTightFigure(angleTime,'Figures/angleTime_modified.pdf')