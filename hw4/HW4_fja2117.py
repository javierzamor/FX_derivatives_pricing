import scipy as sp
import pandas as pd
def main():
    T = pd.read_csv('fx_vol_data.csv')
    label=['1w','1m','6m','1y']
    print(' ------- ---------------------------- ------------------------------- ')
    print('|       |     StdDev of intraday     | StdDv of intraday "perceived" |')
    print('|       |    volatility variation    |   volatility variation of a   |')
    print('| Tenor | (volatility of volatility) |  Vega Neutral Hedged AUDJPY   |')
    print('|       |  of AUDJPY currency pair   |       options portfolio       |')
    print('|-------|----------------------------|-------------------------------|')
    for i in range(4):
        delta_vol=sp.diff(T.iloc[:,9+i],n=1,axis=0)
        crr = sp.array((T.iloc[:,1+i]**2+T.iloc[:,5+i]**2-T.iloc[:,9+i]**2)/(2*T.iloc[:,1+i]*T.iloc[:,5+i]))
        vv=pd.Series(sp.concatenate((crr[0:1],crr[0:-1]),axis=0))
        err_est_cor=sp.sqrt(T.iloc[:,1+i]**2+T.iloc[:,5+i]**2-2*vv*T.iloc[:,1+i]*T.iloc[:,5+i])-T.iloc[:,9+i]
        print('|  '+label[i]+'\t|\t\t'+str(delta_vol.std()) + '\t\t |\t\t\t' + str(err_est_cor.std())+'\t\t |')
    print(' ------- ---------------------------- ------------------------------- ')

if __name__ == '__main__':
    main()
