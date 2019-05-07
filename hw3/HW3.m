function HW3()
    clear;
    clc;
    close all;
    fclose all;
    
    MM = csvread('results.csv');
       
    spot = 1;
    vol_atm  = 0.08;
    dRR_25 = 0.01;
    dRR_10 = 0.018;
    dBF_25 = 0.0025;
    dBF_10 = 0.0080;
    T = 0.5;

    vol=zeros(1,5);
    vol(1)=vol_atm+dBF_10-0.5*dRR_10 ;
    vol(2)=vol_atm+dBF_25-0.5*dRR_25 ;
    vol(3)=vol_atm;
    vol(4)=vol_atm+dBF_25+0.5*dRR_25 ;
    vol(5)=vol_atm+dBF_10+0.5*dRR_10 ;
    
    K=zeros(1,5);
    K(1)  = spot * exp(vol(1) * vol(1) *0.5* T + vol(1) * sqrt(T) * norminv(0.10)); %10p
    K(2)  = spot * exp(vol(2) * vol(2) * 0.5*T + vol(2) * sqrt(T) * norminv(0.25)); %25p
    K(3) = spot * exp(vol(3) * vol(3) * 0.5*T);
    K(4)  = spot * exp(vol(4) * vol(4) * 0.5*T - vol(4) * sqrt(T) * norminv(0.25)); %25c
    K(5)  = spot * exp(vol(5) * vol(5) * 0.5*T - vol(5) * sqrt(T) * norminv(0.10)); %10c
        
    F=[0.01, 10];
    K_min=K(1)*exp(-F*vol(1)*sqrt(T));
    K_max=K(5)*exp(F*vol(5)*sqrt(T));
    
    [K1,v1] = fit_vol(vol,K,K_min(1),K_max(1));
    [K2,v2] = fit_vol(vol,K,K_min(2),K_max(2));
       
    
    disp('');
    h=figure();
    set(h,'windowstyle','docked');
    subplot(1,2,1);
    plot(K,vol,'rs');
    hold on;
    
    plot(MM(:,1),MM(:,2),'k');
    plot(K1,v1,'go');
    grid on;
    hold off;
    
    subplot(1,2,2);
    plot(K,vol,'rs');
    hold on;
    plot(MM(:,3),MM(:,4),'k');
    plot(K2,v2,'go');
    grid on;
    hold off;
    disp('');

end




function [strk,vol]=fit_vol(Y,K,K_min,K_max)

    X = [K_min K K_max];
    X2=X.^2;
    X3=X.^3;
   
    M=[-1	X3(1)	X2(1)	X(1)	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	X3(2)	X2(2)	X(2)	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	X3(2)	X2(2)	X(2)	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	X3(3)	X2(3)	X(3)	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	X3(3)	X2(3)	X(3)	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	X3(4)	X2(4)	X(4)	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	X3(4)	X2(4)	X(4)	1	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	X3(5)	X2(5)	X(5)	1	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	X3(5)	X2(5)	X(5)	1	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	X3(6)	X2(6)	X(6)	1	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	X3(6)	X2(6)	X(6)	1	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	X3(7)	X2(7)	X(7)	1	-1;
        0	-3*X2(1) -2*X(1) -1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	3*X2(2)	2*X(2)	1	0	-3*X2(2)	-2*X(2)	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	3*X2(3)	2*X(3)	1	0	-3*X2(3)	-2*X(3)	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	3*X2(4)	2*X(4)	1	0	-3*X2(4)	-2*X(4)	-1	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	3*X2(5)	2*X(5)	1	0	-3*X2(5)	-2*X(5)	-1	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3*X2(6)	2*X(6)	1	0	-3*X2(6)	-2*X(6)	-1	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3*X2(7)	2*X(7)	1	0	0;
        0	-6*X(1)	-2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	6*X(2)	2	0	0	-6*X(2)	-2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	6*X(3)	2	0	0	-6*X(3)	-2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	6*X(4)	2	0	0	-6*X(4)	-2	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	6*X(5)	2	0	0	-6*X(5)	-2	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	6*X(6)	2	0	0	-6*X(6)	-2	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	6*X(7)	2	0	0	0];

    n=[0 Y(1) Y(1) Y(2) Y(2) Y(3) Y(3) Y(4) Y(4) Y(5) Y(5) 0	0	0	0	0	0	0	0	0	0	0	0	0	0 0]';
    
    coef=M\n;
    coef=[0;0;0;coef(1:end-1);0;0;0;coef(end)];
    coef=reshape(coef,4,[]);
    a = coef(1,:);
    b = coef(2,:);
    c = coef(3,:);
    d = coef(4,:);

    max_strk=K(end)+0.50*(K(end)-K(1));
    min_strk=K(1)-0.25*(K(end)-K(1));
    num_points=40;
    dstrike=(max_strk - min_strk) / (num_points - 1);
    vol=zeros(num_points,1);
    strk=zeros(num_points,1);
    for i=0:(num_points-1)
        strk(i+1)=min_strk + i * dstrike;
        ind=max(find(X<=strk(i+1)));
        if isempty(ind)
            ind=0;
        end
 
        vol(i+1)=a(ind+1)*(strk(i+1))^3+ b(ind+1)*(strk(i+1))^2+c(ind+1)*strk(i+1)+d(ind+1) ;
    end
    disp('');
    
end
