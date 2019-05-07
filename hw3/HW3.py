import scipy
import scipy.stats
import matplotlib.pyplot as plt

def fit_vol(Y,K,F,T):
    K_min=K[0]*scipy.exp(-F*Y[0]*scipy.sqrt(T))
    K_max=K[4]*scipy.exp(F*Y[4]*scipy.sqrt(T))
    X=scipy.concatenate((K_min,K,K_max),axis=0)
    X2=scipy.multiply(X,X)
    X3=scipy.multiply(X2,X)

    M=scipy.matrix([[-1,float(X3[0]),float(X2[0]),float(X[0]),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,	float(X3[1]),float(X2[1]),float(X[1]),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,float(X3[1]),float(X2[1]),float(X[1]),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,float(X3[2]),float(X2[2]),float(X[2]),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,float(X3[2]),float(X2[2]),float(X[2]),1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,float(X3[3]),float(X2[3]),float(X[3]),1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[3]),float(X2[3]),float(X[3]),1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[4]),float(X2[4]),float(X[4]),1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[4]),float(X2[4]),float(X[4]),1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[5]),float(X2[5]),float(X[5]),1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[5]),float(X2[5]),float(X[5]),1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,float(X3[6]),float(X2[6]),float(X[6]),1,-1],
        [0,-3*float(X2[0]),-2*float(X[0]),-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,3*float(X2[1]),2*float(X[1]),1,0,-3*float(X2[1]),-2*float(X[1]),-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,3*float(X2[2]),2*float(X[2]),1,0,-3*float(X2[2]),-2*float(X[2]),-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,3*float(X2[3]),2*float(X[3]),1,0,-3*float(X2[3]),-2*float(X[3]),-1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,3*float(X2[4]),2*float(X[4]),1,0,-3*float(X2[4]),-2*float(X[4]),-1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3*float(X2[5]),2*float(X[5]),1,0,-3*float(X2[5]),-2*float(X[5]),-1,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3*float(X2[6]),2*float(X[6]),1,0,0],
        [0,-6*float(X[0]),-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,6*float(X[1]),2,0,0,-6*float(X[1]),-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,6*float(X[2]),2,0,0,-6*float(X[2]),-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,6*float(X[3]),2,0,0,-6*float(X[3]),-2,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,6*float(X[4]),2,0,0,-6*float(X[4]),-2,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6*float(X[5]),2,0,0,-6*float(X[5]),-2,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6*float(X[6]),2,0,0,0]])

    n=scipy.matrix([0,float(Y[0]),float(Y[0]),float(Y[1]),float(Y[1]),float(Y[2]),float(Y[2]),float(Y[3]),float(Y[3]),float(Y[4]),float(Y[4]),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    n=scipy.matrix.transpose(n)
    coef=M.I*n

    z3=scipy.matrix(scipy.zeros((3,1)))
    coef=scipy.concatenate((z3,coef[0:-1],z3,coef[-1]),axis=0)
    coef=scipy.reshape(coef,(8,4))
    a=coef[:,0]
    b=coef[:,1]
    c=coef[:,2]
    d=coef[:,3]

    max_strk=float(K[-1]+0.50*(K[-1]-K[0]))
    min_strk=float(K[0]-0.25*(K[-1]-K[0]))
    num_points=40
    step=(max_strk-min_strk)/(num_points-1)
    vol=scipy.matrix(scipy.zeros((num_points,1)))
    strk=scipy.arange(min_strk,max_strk+step,step)
    for cont in range(scipy.size(strk)):
        inx=(X<=strk[cont])
        if scipy.all(inx==0):
            inx=0
        else:
            inx=inx.nonzero()
            inx=max(inx[0])+1
        strkk=float(strk[cont])
        vol[cont]=float(a[inx])*(strkk**3)+float(b[inx])*(strkk**2)+float(c[inx])*strkk+float(d[inx])

    return strk,vol

def main():

    spot=1
    vol_atm=0.08
    dRR_25=0.01
    dRR_10=0.018
    dBF_25=0.0025
    dBF_10=0.0080
    T=0.5
    F=[0.01, 10]

    vol=scipy.matrix(scipy.zeros((5,1)))
    vol[0]=vol_atm+dBF_10-0.5*dRR_10
    vol[1]=vol_atm+dBF_25-0.5*dRR_25
    vol[2]=vol_atm
    vol[3]=vol_atm+dBF_25+0.5*dRR_25
    vol[4]=vol_atm+dBF_10+0.5*dRR_10

    K=scipy.matrix(scipy.zeros((5,1)))
    K[0]=spot*scipy.exp(vol[0]*vol[0]*0.5*T+vol[0]*scipy.sqrt(T)*scipy.stats.norm.ppf(0.10))
    K[1]=spot*scipy.exp(vol[1]*vol[1]*0.5*T+vol[1]*scipy.sqrt(T)* scipy.stats.norm.ppf(0.25))
    K[2]=spot*scipy.exp(vol[2]*vol[2]*0.5*T)
    K[3]=spot*scipy.exp(vol[3]*vol[3]*0.5*T-vol[3]*scipy.sqrt(T)*scipy.stats.norm.ppf(0.25))
    K[4]=spot*scipy.exp(vol[4]*vol[4]*0.5*T-vol[4]*scipy.sqrt(T)*scipy.stats.norm.ppf(0.10))

    [K1,v1] = fit_vol(vol,K,F[0],T)
    [K2,v2] = fit_vol(vol,K,F[1],T)

    plt.plot(K1, v1)
    plt.plot(K2, v2)
    plt.plot(K,vol,'go')
    plt.xlabel('Strike')
    plt.ylabel('Volatility')
    plt.grid()
    plt.show(block=False)

if __name__ == '__main__':
     main()


