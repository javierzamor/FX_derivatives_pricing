function HW4()
    
    clc;
    clear;
    close all;
    fclose all;
    
    T=readtable('fx_vol_data.csv');
    
    TT=[T.(10) T.(2) T.(6) T.(11) T.(3) T.(7) T.(12) T.(4) T.(8) T.(13) T.(5) T.(9)];
    
    label=['1w';'1m';'6m';'1y'];
    
    fprintf(' ------- ---------------------------- --------------------- --------- \n'); 
    fprintf('|       |     StdDev of intraday     | StdDv of intraday "perceived" |\n');
    fprintf('|       |    volatility variation    |   volatility variation of a   |\n');
    fprintf('| Tenor | (volatility of volatility) |  Vega Neutral Hedged AUDJPY   |\n');
    fprintf('|       |  of AUDJPY currency pair   |       options portfolio       |\n')
    fprintf('|-------|----------------------------|-------------------------------|\n');
    for i=0:3
        c0=TT(:,i*3+1);
        c1=TT(:,i*3+2);
        c2=TT(:,i*3+3);
        
        delta_vol=diff(c0);
        std_vol_vol=nanstd(delta_vol);
        
        crr = (c1.^2+c2.^2-c0.^2)./(2.*c1.*c2);
        v_tmp=[crr(1);crr(1:end-1)];
        est_c0=sqrt(c1.^2+c2.^2-2.*v_tmp.*c1.*c2);
        err_est_cor=est_c0-c0;
        std_err_est_cor=nanstd(err_est_cor);
        
        
        fprintf('|  %s   |\t\t\t%5.5f\t\t\t |\t\t\t\t%5.5f\t\t\t |\n',label(i+1,:),std_vol_vol,std_err_est_cor);
    end
    fprintf(' ------- ---------------------------- ------------------------------- \n');  
end

