
import scipy

def main():
    vol1=0.01
    vol2=0.008
    b1=0.5
    b2=0.1
    rho=-0.4
    T1=0.25
    T2=1
    Q=0.03
    dt=1e-3
    numruns=100000
    scipy.random.seed(1)

    rr=scipy.zeros((4,6))
    rr[0,]=[0.1, 0.25, 0.5, 0.75, 1, 2]
    for contT in range(0,6):
        T=rr[0,contT]
        for hedge in range(0,3):
            dz1=scipy.random.normal(0,scipy.sqrt(dt), numruns)
            dz2=rho*dz1+scipy.sqrt(1-rho**2)*scipy.random.normal(0,scipy.sqrt(dt), numruns)
            dQT=vol1*scipy.exp(-b1*T)*dz1 +vol2*scipy.exp(-b2*T)*dz2
            dQ1=vol1*scipy.exp(-b1*T1)*dz1 +vol2*scipy.exp(-b2*T1)*dz2
            dQ2=vol1*scipy.exp(-b1*T2)*dz1 +vol2*scipy.exp(-b2*T2)*dz2
            pnl=scipy.exp(-(Q+dQT)*T)-scipy.exp(-Q*T)
            pnl=pnl-N1(vol1,vol2,b1,b2,T,T1,T2,Q,hedge)*(scipy.exp(-(Q+dQ1)*T1)-scipy.exp(-Q*T1))
            pnl=pnl-N2(vol1,vol2,b1,b2,T,T1,T2,Q,hedge)*(scipy.exp(-(Q+dQ2)*T2)-scipy.exp(-Q*T2))
            rr[hedge+1,contT]=pnl.std()
    scipy.set_printoptions(precision=3)
    print(rr)

def N1(vol1,vol2,b1,b2,T,T1,T2,Q,hedge):
    Not1=0
    dT=T2-T1
    if hedge == 1:
        if T<=T1:
            Not1=T*scipy.exp(-Q*(T-T1))/T1
        elif T>T1 and T<T2:
            Not1=(T2-T)*T*scipy.exp(-Q*(T-T1))/(dT*T1)
    elif hedge==2:
        dz1=-scipy.exp(b1*T2-b2*dT)/(vol1*(1-scipy.exp((b1-b2)*dT)))
        dz2=scipy.exp(b2*T1)/(vol2*(1-scipy.exp((b1-b2)*dT)))
        dQT=vol1*scipy.exp(-b1*T)*dz1+vol2*scipy.exp(-b2*T)*dz2
        Not1=dQT*T/T1*scipy.exp(-Q*(T-T1))
    return Not1

def N2(vol1,vol2,b1,b2,T,T1,T2,Q,hedge):
    Not2=0
    dT=T2-T1
    if hedge == 1:
        if T>=T2:
            Not2=T*scipy.exp(-Q*(T-T2))/T2
        elif T>T1 and T<T2:
            Not2=(T-T1)*T*scipy.exp(Q*(T2-T))/(dT*T2)
    elif hedge==2:
        dz1=-scipy.exp(b1*T1+b2*dT)/(vol1*(1-scipy.exp(-(b1-b2)*dT)))
        dz2=scipy.exp(b2*T2)/(vol2*(1-scipy.exp(-(b1-b2)*dT)))
        dQT=vol1*scipy.exp(-b1*T)*dz1+vol2*scipy.exp(-b2*T)*dz2
        Not2=dQT*T*scipy.exp(Q*(T2-T))/T2
    return Not2

if __name__=="__main__":
    main()


