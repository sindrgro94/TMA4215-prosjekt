function [tEvent,totEvent] = hermiteInterpolationPecker(t,y,finding)
    h = t(2)-t(1);
    P0 = y(1,1);
    P3 = y(1,2);
    P1 = P0+h*y(2,1)/3;
    P2 = P3-h*y(2,2)/3;
    BO = @(theta) P0*(1-theta)^3+3*P1*(1-theta)^2*theta+...
        3*P2*(1-theta)*theta^2+P3*theta^3;
    BOder = @(theta) -3*P0*(1-theta)^2+3*P3*theta^2+...
        3*P1*((1-theta)^2-2*theta*(1-theta))+3*P2*(2*(1-theta)*theta-theta^2);
    

    
    %Newtons method
    thetaN = 0.5;
    while (abs(BO(thetaN)-finding))>10^-12
         thetaN = thetaN - ((BO(thetaN)-finding)/BOder(thetaN));
    end
    tEvent = t(1)+thetaN*h;
    yEvent = [BO(thetaN); BOder(thetaN)/h];
    
    Q0 = y(3,1);
    Q3 = y(3,2);
    Q1 = Q0+h*y(4,1)/3;
    Q2 = Q3-h*y(4,2)/3;
    Bz = @(theta) Q0*(1-theta)^3+3*Q1*(1-theta)^2*theta+...
        3*Q2*(1-theta)*theta^2+Q3*theta^3;
    Bzder = @(theta) -3*Q0*(1-theta)^2+3*Q3*theta^2+...
        3*Q1*((1-theta)^2-2*theta*(1-theta))+3*Q2*(2*(1-theta)*theta-theta^2);
    zEvent = [Bz(tEvent);Bzder(tEvent)];
    
    totEvent = [yEvent;zEvent];
    
    
    
    
    
    
    
    
end