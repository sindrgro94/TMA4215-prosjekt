function [tEvent,yEvent] = hermiteInterpolation(t,y)
    h = t(2)-t(1);
    P0 = y(1,1);
    P3 = y(1,2);
    P1 = P0+h*y(2,1)/3;
    P2 = P3-h*y(2,2)/3;
    B = @(theta) P0*(1-theta)^3+3*P1*(1-theta)^2*theta+...
        3*P2*(1-theta)*theta^2+P3*theta^3;
    Bder = @(theta) -3*P0*(1-theta)^2+3*P3*theta^2+...
        3*P1*((1-theta)^2-2*theta*(1-theta))+3*P2*(2*(1-theta)*theta-theta^2);
    %Newtons method
    thetaN = 0.5;
    while abs(B(thetaN))>10^-12
         thetaN = thetaN - (B(thetaN)/Bder(thetaN));
    end
    tEvent = t(1)+thetaN*h;
    yEvent = [B(thetaN); Bder(thetaN)/h];
end