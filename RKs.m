function [t, y, iflag, nfun, njac] = RKs(f, jac, t0, tend, y0, Tol, h0,eventLocator)
    nfun = 0;
    njac = 0;
    steps = (tend-t0)/h0*5; %estimation of steps
    t = zeros(1,steps);
    y = zeros(length(y0),steps);
    h = h0;
    y(:,1) = y0;
    t(1) = t0;
    step = 1;
    newtonTol = 0.1*Tol;
    P = 0.8; %pessimist factor
    stepTol = 10^(-15); %no idea what this should be.
    while t(step)<tend
        [tnext, ynext, le, iflag, funEv, jacEv] = onestep(f,jac,t(step),y(:,step),h,newtonTol);
        nfun = nfun + funEv;
        njac = njac + jacEv;
        %[tnext, ynext, le, iflag] = onestep_old(f,jac,t(step),y(:,step),h,newtonTol);
        if iflag == 1 && le<Tol
            t(step+1) = tnext;
            y(:,step+1) = ynext;
            h = P*(Tol/le)^(1/3)*h;
            step = step + 1;
            if t(step)+h>tend
                h = tend-t(step);
            end
            if eventLocator{1} == true
                if ynext >= eventLocator{2}
                    return
                end
            end
        elseif h<stepTol
                iflag = -1;
                return
        else
            h = 0.5*h;
        end
end