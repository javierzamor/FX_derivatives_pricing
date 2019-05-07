function HW_2()
    clear
    format shortg;
    sig1=0.01;
    sig2=0.008;
    b1=0.5;
    b2=0.1;
    rho=-0.4;
    T1=0.25;
    T2=1;
    Q=0.03;
    dt=1e-3;
    num_runs=1e6;
    rng(1);

    rr=zeros(4,6);
    rr(1,:)=[0.1,0.25,0.5,0.75,1,2];

    for contT=1:6
        for model=0:2
            T=rr(1,contT);
            dz1=random('Normal',0,sqrt(dt),[num_runs,1]);
            dz2=rho*dz1+sqrt(1-rho^2)*random('Normal',0,sqrt(dt),[num_runs,1]);
            dQT = sig1 * exp(-b1 * T)  * dz1 + sig2 * exp(-b2 * T)  * dz2;
            dQ1 = sig1 * exp(-b1 * T1) * dz1 + sig2 * exp(-b2 * T1) * dz2;
            dQ2 = sig1 * exp(-b1 * T2) * dz1 + sig2 * exp(-b2 * T2) * dz2;
            pnls= (exp(-(Q + dQT) * T)  - exp(-Q * T));
            pnls=pnls-N1(sig1,sig2,b1,b2,T,T1,T2,Q,model) * (exp(-(Q + dQ1) * T1) - exp(-Q * T1));
            pnls=pnls-N2(sig1,sig2,b1,b2,T,T1,T2,Q,model) * (exp(-(Q + dQ2) * T2) - exp(-Q * T2));
            rr(model+2,contT)=std(pnls)*1e4;
        end
        disp('');
    end
    
    rr
    disp('');
end

function Not1=N1(sig1,sig2,b1,b2,T,T1,T2,Q,model)
    dT=T2-T1;
    Not1=0;
    if model ==1
        if T<=T1
            Not1=T*exp(-Q*(T-T1))/T1;
        elseif T>T1 && T<T2
            Not1=T*(T2-T)*exp(-Q*(T-T1))/(T1*dT);
        end
    elseif model==2
        dz1 = -exp(b1*T2-b2*dT)/(sig1*(1-exp((b1-b2)*dT)));
        dz2 =  exp(b2*T1)/(sig2*(1-exp((b1-b2)*dT)));
        dQT=sig1*exp(-b1*T)*dz1+sig2*exp(-b2*T)*dz2;
        Not1 = dQT*T*exp(-Q*(T-T1))/T1 ;
    end


end

function Not2=N2(sig1,sig2,b1,b2,T,T1,T2,Q,model)
    dT=T2-T1;
    Not2=0;
    if model==1
        if T>=T2
            Not2=T*exp(-Q*(T-T2))/T2;
        elseif T>T1 && T<T2
            Not2=T*(T-T1)*exp(-Q*(T-T2))/(T2*dT) ;
        end
    elseif model==2
        dz1=-exp(b1*T1+b2*dT)/(sig1*(1-exp(-(b1-b2)*dT)));
        dz2= exp(b2*T2)/(sig2*(1-exp(-(b1-b2)*dT)));
        dQT=sig1*exp(-b1*T)*dz1+sig2*exp(-b2*T)*dz2;
        Not2=dQT*T*exp(-Q*(T-T2))/T2;
    end
end
